"""util_infer.py

"""
from __future__ import absolute_import
from __future__ import division

import numpy
import multiprocessing
import itertools
import datetime
import os
import math
try:
    import cPickle as pickle
except:
    import pickle

import cmpy
from cmpy.math import log, logaddexp
import cmpy.inference.bayesianem as bayesem
import cmpy.orderlygen.pyicdfa as pyidcdfa

from .util_general import deltatime_format
from .util_general import read_evidence_file
from .util_general import read_machines_file
from .util_general import read_probabilities_file
from .util_general import read_sample_dir
from .util_general import write_evidence_file
from .util_general import write_probabilities_file

def calculate_probabilities(dbdir, inferemdir, beta, penalty):
    """A function to calculate the probabilities for a file containing
    machine/model evidence values.

    Parameters
    ----------
    dbdir : str
        Base directory for the database (directory).
    inferemdir : str
        Sub-directory for log_evidence file.
    beta : float
        Strength of exponential penalty.
    penalty : str
        Penalty can be for: num_states or num_edges.
    
    Returns
    -------
    (summary_str, model_probabilities) : (str, dict)
        A summary of function operations and a dictionary containing model
        probabilities.

    """
    summary = []

    # location of inferEM instances
    inferdir = os.path.join(dbdir, inferemdir) 

    # add to summary list
    summary.append("* using DB: {}\n".format(dbdir))
    summary.append("  - subdir: {}\n\n".format(inferdir))
    
    # start processing...
    script_start = datetime.datetime.now()
    script_start_str = script_start.strftime("%H:%M:%S %D")
    summary.append(" -- start time   : {}\n".format(script_start_str))

    # open evidence dictionary from disk
    evidence_dict = read_evidence_file(inferdir)
    
    # process log evidence terms
    log_evidence_total = log(0)
    log_evidence = {}
    for mname, minfo in evidence_dict.items():
        # find penalty type and amount
        if penalty == 'num_states':
            penalty_term = minfo['nodes']
        elif penalty == 'num_edges':
            penalty_term = minfo['edges']

        log_evidence[mname] = minfo['log_evidence'] - beta*penalty_term
        log_evidence_total = logaddexp(log_evidence_total,log_evidence[mname])
    
    # calculate probabilities
    model_probabilities = {}
    for em in log_evidence.keys():
        # divide log_evidence by normalization
        temp_log_prob = -log_evidence_total + log_evidence[em]
        
        # calculate probablirt
        model_probabilities[em] = math.exp(temp_log_prob)
    
    # check for existing output file before doing work
    ftemp = "probabilities_beta-{:f}_penalty-{:s}".format(beta, penalty)
    fname = os.path.join(inferdir, ftemp)
    if os.path.exists(fname):
        raise Exception("\nOutput file {} exists!\n".format(fname))

    # save probabilities to file
    write_probabilities_file(model_probabilities, fname)
    write_probabilities_file(model_probabilities, fname + '_top100', 100)

    # end processing...
    script_end = datetime.datetime.now()
    script_end_str = script_end.strftime("%H:%M:%S %D")
    summary.append(" -- end time     : {}\n".format(script_end_str))
    time_diff = script_end-script_start
    time_diff_str = deltatime_format(time_diff)
    summary.append(" -- compute time : {}\n\n".format(time_diff_str))

    summary_str = ''.join(summary)

    return (summary_str, model_probabilities)

def create_model_string(model):
    """Create a string representation of model for writing to file."""
    import cmpy
    from cmpy.machines import MealyHMM
    
    ## default assumption type is uHMM
    model_type = 'uHMM'

    # need to switch to MealyHMM class for use of is_minimal_unifilar() method
    model.__class__ = MealyHMM  
    if model.is_minimal_unifilar():
        # topological eM
        model_type = 'topEM'

    model_topology = model.to_string().replace('\n','')

    # create string
    model_str = "{},{},{},{},{}\n".format(str(model),
                                          model_type,
                                          len(model.nodes()),
                                          len(model.edges()),
                                          model_topology)
    
    return model_str

def create_machine_file(filename, A, N_list, em_min, nmax, nprocs):
    """Create a file with candidate model information.

    Parameters
    ----------
    filename : str
        Complete path for the machine file.
    A : int
        Alphabet size.
    N_list : list (ints)
        A list of node (model) sizes to consider.
    em_min : str
        Consider just topological eMs or all uHMM topologies?
    nmax : int
        Maximum number of machine to load into RAM before processing.
    nprocs : int
        Number of simultaneous processes for mutiprocessing.
        
    Returns
    -------
    summary_str : str
        A summary with compute time information.
    
    """
    # start processing...
    script_start = datetime.datetime.now()
    script_start_str = script_start.strftime("%H:%M:%S %D")
    summary = []
    summary.append(" -- start time   : {}\n".format(script_start_str))

    # open file and write header information
    f = open(filename, 'w')
    f.write("# candidate machines\n")
    f.write("# - alphabet size: {}\n".format(A))
    f.write("# - number of states: {}\n".format(N_list))
    if em_min == 'none':
        topologies = 'all uHMMS'
    else:
        topologies = 'only topological eMs'
    f.write("# - topologies: {}\n".format(topologies))

    # column header
    out_str = "{},{},{},{},{}\n".format('name', 'HMM type',
                                        'states', 'edges', 'topology')
    f.write(out_str)
    
    # initialize model topology iterator
    model_iter = bayesem.LibraryGenerator(A, N_list, em_min)

    # create list of models
    model_counter = 0
    model_strings = {}
    models = []
    for m in model_iter:
        models.append(m)
        model_counter += 1

        if len(models) == nmax:
            # farm out models enumerated thus far
            model_strings.update(mp_model_strings(models, nprocs))
            models = []
        else:
            pass

    # need to check for left over machines
    model_strings.update(mp_model_strings(models, nprocs))

    for m_name in sorted(model_strings.iterkeys()):
        f.write(model_strings[m_name])
    
    f.close()

    # end processing...
    script_end = datetime.datetime.now()
    script_end_str = script_end.strftime("%H:%M:%S %D")
    summary.append(" -- end time     : {}\n".format(script_end_str))
    time_diff = script_end-script_start
    time_diff_str = deltatime_format(time_diff)
    summary.append(" -- compute time : {}\n\n".format(time_diff_str))

    # write out number of topologies
    summary.append("\n** Number of topologies: {}\n\n".format(model_counter))

    summary_str = ''.join(summary)

    return summary_str

def create_machine_posterior_file(dbdir, range, data, nprocs):
    """A fuction to generate a text file consisting of log_evidence for the
    posterior.

    Parameters
    ----------
    dbdir : str
        The database (directory) name
    range : tuple (int, int)
        A tuple of integers provide the data range to consider
    data : list
        The complete data series -- assumes appropriate slice already made
    nprocs : int
        Number of simultaneous processes for mutiprocessing

    Returns
    -------
    summary_str : str
        A summary of the inference done for print to stdout or file.

    """
    cwd = os.getcwd()
    summary = []

    # mkdir inferEM_pt1-pt2 for posterior
    pt1 = range[0]
    pt2 = range[1]
    inferdir = os.path.join(dbdir,"inferEM_{:d}-{:d}".format(pt1, pt2))
    summary.append("* making subdirectory: {}\n".format(inferdir))
    os.mkdir(inferdir)

    # change to inferEM database directory
    os.chdir(inferdir)
    summary.append("* using DB: {} subdir: {}\n\n".format(dbdir, inferdir))

    # start inference...
    script_start = datetime.datetime.now()
    script_start_str = script_start.strftime("%H:%M:%S %D")
    summary.append(" -- start time   : {}\n".format(script_start_str))

    # open machines file
    machines = read_machines_file(os.path.join(cwd,dbdir))

    # do serious processing....
    evidence_dict = mp_machine_evidence(machines, data, nprocs)

    # write the evidence dictionary to a file
    write_evidence_file(evidence_dict, 'log_evidence')
    write_evidence_file(evidence_dict, 'log_evidence_top100', 100)

    # end processing...
    script_end = datetime.datetime.now()
    script_end_str = script_end.strftime("%H:%M:%S %D")
    summary.append(" -- end time     : {}\n".format(script_end_str))
    time_diff = script_end-script_start
    time_diff_str = deltatime_format(time_diff)
    summary.append(" -- compute time : {}\n\n".format(time_diff_str))

    summary_str = ''.join(summary)

    # return to cwd
    os.chdir(cwd)

    return summary_str

def create_machine_prior_file(dbdir, nprocs):
    """A fuction to generate a text file consisting of log_evidence file for
    prior.

    Parameters
    ----------
    dbdir : str
        The database (directory) name.
    nprocs : int
        Number of simultaneous processes for mutiprocessing.

    Returns
    -------
    summary_str : str
        A summary of the inference done for print to stdout or file.

    """
    cwd = os.getcwd()
    summary = []

    # mkdir inferEM_0-0 for prior
    inferdir = os.path.join(dbdir,"inferEM_{:d}-{:d}".format(0, 0))
    summary.append("* making subdirectory: {}\n".format(inferdir))
    os.mkdir(inferdir)
    
    # change to inferEM database directory
    os.chdir(inferdir)
    summary.append("* using DB: {} subdir: {}\n\n".format(dbdir, inferdir))

    # start inference...
    script_start = datetime.datetime.now()
    script_start_str = script_start.strftime("%H:%M:%S %D")
    summary.append(" -- start time   : {}\n".format(script_start_str))

    # open machines file
    machines = read_machines_file(os.path.join(cwd,dbdir)) 

    # do serious processing....
    evidence_dict = mp_machine_evidence(machines, None, nprocs)

    # write the evidence dictionary to a file
    write_evidence_file(evidence_dict, 'log_evidence')
    write_evidence_file(evidence_dict, 'log_evidence_top100', 100)
    
    # end processing...
    script_end = datetime.datetime.now()
    script_end_str = script_end.strftime("%H:%M:%S %D")
    summary.append(" -- end time     : {}\n".format(script_end_str))
    time_diff = script_end-script_start
    time_diff_str = deltatime_format(time_diff)
    summary.append(" -- compute time : {}\n\n".format(time_diff_str))

    summary_str = ''.join(summary)

    # return to cwd
    os.chdir(cwd)

    return summary_str

def create_sample_summary_file(db_dir, sample_dir, nprocs):
    """Create a file with properties of sample machines from prior or
    posterior.

    Parameters
    ----------
    db_dir : str
        Name of the root database directory.
    sample_dir : str
        Name of the subdirectory containing the sample machine in pickled form.
    
    Returns
    -------
    summary : str
        A summary with compute time information.

    """
    # start processing...
    summary = []
    script_start = datetime.datetime.now()
    script_start_str = script_start.strftime("%H:%M:%S %D")
    summary.append(" -- start time   : {}\n".format(script_start_str))

    # get file 
    sample_files = read_sample_dir(db_dir, sample_dir)

    # open file for writing
    foutname = db_dir + '/' + sample_dir + '.dat'
    fout = open(foutname, 'w')

    # header
    temp = ["Filename",
            "Machine_ID",
            "Recurrent_eM",
            "State_Entropy",
            "Entropy_Rate",
            "Number_Edges",
            "Number_States\n"
            ]        
    fout.write(','.join(temp))

    # change dir
    startdir = os.getcwd()
    os.chdir(os.path.join(db_dir, sample_dir))

    # send to multiprocessing routine
    sample_info = mp_process_samples(sample_files, nprocs)

    # write output to file
    for sfile in sorted(sample_info.keys()):
        fout.write(','.join(sample_info[sfile]))

    fout.close()

    # change back to startdir
    os.chdir(startdir)
    
    script_end = datetime.datetime.now()
    script_end_str = script_end.strftime("%H:%M:%S %D")
    summary.append(" -- end time     : {}\n".format(script_end_str))
    time_diff = script_end-script_start
    time_diff_str = deltatime_format(time_diff)
    summary.append(" -- compute time : {}\n\n".format(time_diff_str))

    summary_str = ''.join(summary)

    return summary_str

def create_sample_machine_file(sample_num, sampler, sampledir, data):
    """A map function for mutiprocessing to consider sampling of machines.
    
    Parameters
    ----------
    sample_num : int
        The sample number -- used for naming output files uniquely.
    sampler : DictionarySampler
        An instance of the DictionarySampler class.
    sampledir : str
        Path to output directory for sampled machines.

    """
    import cmpy
    import cmpy.orderlygen.pyicdfa as pyidcdfa
    import cmpy.inference.bayesianem as bayesem

    # sample a model name
    sample_mname = sampler()
    
    # create filename for sampled model name
    sample_file = ''.join(['inferEM_', sample_mname, '.pickle'])

    # - number of states
    n = int(sample_mname.split('_')[0][1:])
    # - alphabet size
    k =  int(sample_mname.split('_')[1][1:])
    # - id
    id =  int(sample_mname.split('_')[2][2:])

    # get machine topology
    machine = pyidcdfa.int_to_machine(id, n, k)
    # set name, with proper ID - n states, k symbols
    mname = "n{}_k{}_id{}".format(n, k, id)
    machine.set_name(mname)

    # generate inferEM instance
    inferem_instance = bayesem.InferEM(machine, data)

    # sample machine (also returns start node, not needed)
    _, em_sample = inferem_instance.generate_sample()
    
    # create file name
    fname = "sample{:05d}.pickle".format(sample_num)
    
    # write to file
    f = open(os.path.join(sampledir,fname), 'w')
    pickle.dump(em_sample, f)
    f.close()

    return sample_num

def machine_log_evidence(machine, data):
    """A function for getting the log evidence for a single model topology and
    data set.  If None is passed for data, get values for prior.

    Paremeters
    ----------
    machine : MealyHMM or RecurrentEpsilonMachine
        Model topolgy to consider
    data : list, None
        Data series to consider for inference
    
    Return
    ------
    (mname, mdict) : tuple
        Tuple with machine name `mname` and dictionary `mdict` containing
        `log_evidence` for the data series provided as well as the number of
        `nodes` and `edges` for the machine.
    
    """
    # pass to InferEM
    infer_temp = bayesem.InferEM(machine,data)

    # extract machine name and other information
    model_info = (str(machine),{'log_evidence': infer_temp.log_evidence(),
                                'nodes': len(machine.nodes()),
                                'edges': len(machine.edges())})

    return model_info

def mp_machine_evidence(machines, data, nprocs):
    """Multiprocessing control for create prior log evidence file."""

    def mp_worker(machines_subset, data, out_q):
        """Worker for mp_machine_evidence."""
        import cmpy
        outdict = {}
        
        for m in machines_subset:
            # split out elements of list
            # em_name, em_type, em_num_states, em_num_edges, em_str
            em = cmpy.machines.from_string(m[4])
            em.set_name(m[0])

            # get log evidence
            mname, minfo = machine_log_evidence(em, data)
            outdict[mname] = minfo
    
        out_q.put(outdict)
    
    # Each process will get 'chunksize' machines, data and out_q for storing
    # output
    out_q = multiprocessing.Queue()
    chunksize = int(math.ceil(len(machines) / float(nprocs)))
    procs = []

    for i in range(nprocs):
        p = multiprocessing.Process(
                target=mp_worker,
                args=(machines[chunksize * i:chunksize * (i + 1)],
                      data,
                      out_q))

        procs.append(p)
        p.start()
    
    # Collect all results into a single result dict
    resultdict = {}
    for i in range(nprocs):
        resultdict.update(out_q.get())

    # Wait for all worker processes to finish
    for p in procs:
        p.join()

    return resultdict

def mp_model_strings(models, nprocs):
    """Multiprocessing control for create machines file."""

    def mp_worker(models, out_q):
        """Worker for mp_model_strings."""
        outdict = {}
        
        for m in models:
            outdict[str(m)] = create_model_string(m)
    
        out_q.put(outdict)
    
    # Each process will get 'chunksize' nums and a queue to put his out
    # dict into
    out_q = multiprocessing.Queue()
    chunksize = int(math.ceil(len(models) / float(nprocs)))
    procs = []

    for i in range(nprocs):
        p = multiprocessing.Process(
                target=mp_worker,
                args=(models[chunksize * i:chunksize * (i + 1)],
                           out_q))

        procs.append(p)
        p.start()
    
    # Collect all results into a single result dict. We know how many dicts
    # with results to expect.
    resultdict = {}
    for i in range(nprocs):
        resultdict.update(out_q.get())

    # Wait for all worker processes to finish
    for p in procs:
        p.join()

    return resultdict

def mp_process_samples(sample_files, nprocs):
    """Multiprocessing control for analyzing directory of sample machines."""

    def mp_worker(focus_sample_files, out_q):
        """Load sample machine and record various properties."""
        import cmpy
        outdict = {}

        ## cycle through files
        for f in focus_sample_files:
            # load pickled machine
            f_em = open(f, 'r')
            em = pickle.load(f_em)
            f_em.close()
            
            # state entropy, Cmu if epsilon-machine
            try:
                cmu_est = em.statistical_complexity()
            except:
                cmu_est = em.stationary_entropy()

            # entropy rate
            hmu_est = em.entropy_rate()
            
            # number of nodes, edges
            num_state = len(em.nodes())
            num_edges = len(em.edges())

            # get model type
            recurr_em = 0
            if type(em) == cmpy.machines.RecurrentEpsilonMachine:
                recurr_em = 1

            # get model id/name
            em_id = str(em).split(':')[-1].strip()

            temp = ["{fn:s}".format(fn=f),
                    "{mn:s}".format(mn=em_id),
                    "{id:d}".format(id=recurr_em),
                    "{cmu:.15e}".format(cmu=cmu_est),
                    "{hmu:.15e}".format(hmu=hmu_est),
                    "{ne}".format(ne=num_edges),
                    "{ns}\n".format(ns=num_state)
                ]
        
            outdict[f] = temp
        
        out_q.put(outdict)
    
    # Each process will get 'chunksize' number of sample to do
    out_q = multiprocessing.Queue()
    chunksize = int(math.ceil(len(sample_files) / float(nprocs)))
    procs = []

    for i in range(nprocs):
        p = multiprocessing.Process(
                target=mp_worker,
                args=(sample_files[chunksize * i:chunksize * (i + 1)],
                      out_q))

        procs.append(p)
        p.start()
    
    # Collect all results into a single result dict
    resultdict = {}
    for i in range(nprocs):
        resultdict.update(out_q.get())

    # Wait for all worker processes to finish
    for p in procs:
        p.join()

    return resultdict

def mp_sample_models(sampledir, sample_nums, sampler, data, nprocs):
    """Multiprocessing control for create machines file."""

    def mp_worker(sample_nums, sampler, sampledir, data, out_q):
        """Produce a set of sameple eMs or uHMMs with the desired set of sample
        numbers."""

        import cmpy
        import cmpy.orderlygen.pyicdfa as pyidcdfa
        import cmpy.inference.bayesianem as bayesem
        
        # create a local cache for InferEM instances
        inferem_cache = {}

        # run through allocated sample numbers for this process
        outdict = {}
        for sn in sample_nums:
            # sample a model name
            sample_mname = sampler()
            
            if inferem_cache.has_key(sample_mname):
                # already processes, use existing InfereM instance
                inferem_instance = inferem_cache[sample_mname]
            else:
                # create new InferEM instance
                # - number of states
                n = int(sample_mname.split('_')[0][1:])
                # - alphabet size
                k =  int(sample_mname.split('_')[1][1:])
                # - id
                id =  int(sample_mname.split('_')[2][2:])

                # get machine topology
                machine = pyidcdfa.int_to_machine(id, n, k)
                # set name, with proper ID - n states, k symbols
                mname = "n{}_k{}_id{}".format(n, k, id)
                machine.set_name(mname)

                # generate inferEM instance
                inferem_instance = bayesem.InferEM(machine, data)

                # cache for later use
                inferem_cache[sample_mname] = inferem_instance

            # sample machine (also returns start node)
            startnode, em_sample = inferem_instance.generate_sample()
            
            # state entropy, Cmu if epsilon-machine
            try:
                cmu_est = em_sample.statistical_complexity()
            except:
                cmu_est = em_sample.stationary_entropy()

            # entropy rate
            hmu_est = em_sample.entropy_rate()
            
            # number of nodes, edges
            num_state = len(em_sample.nodes())
            num_edges = len(em_sample.edges())

            # get model type
            recurr_em = 'uHMM'
            if type(em_sample) == cmpy.machines.RecurrentEpsilonMachine:
                recurr_em = 'topEM'

            # get model id/name
            em_id = str(em_sample)

            # get str form for sample
            em_str = em_sample.to_string().replace('\n','')

            temp = ["{:d}".format(sn),
                    "{:s}".format(sample_mname),
                    "{:s}".format(recurr_em),
                    "{cmu:.15e}".format(cmu=cmu_est),
                    "{hmu:.15e}".format(hmu=hmu_est),
                    "{ne}".format(ne=num_edges),
                    "{ns}".format(ns=num_state),
                    "{}".format(startnode),
                    "{top}\n".format(top=em_str)
                ]
        
            outdict[sn] = temp
        
        out_q.put(outdict)

    # Each process will get 'chunksize' number of samples to do
    out_q = multiprocessing.Queue()
    chunksize = int(math.ceil(len(sample_nums) / float(nprocs)))
    procs = []

    for i in range(nprocs):
        p = multiprocessing.Process(
                target=mp_worker,
                args=(sample_nums[chunksize * i:chunksize * (i + 1)],
                      sampler,
                      sampledir,
                      data,
                      out_q))

        procs.append(p)
        p.start()
    
    # Collect all results into a single result dict. We know how many dicts
    # with results to expect.
    resultdict = {}
    for i in range(nprocs):
        resultdict.update(out_q.get())
    
    # Wait for all worker processes to finish
    for p in procs:
        p.join()

    return resultdict

def sample_machines(dbdir, inferemdir, modelprobs, num_sample, data, nprocs):
    """A function to generate machine samples from the prior or posterior.

    Parameters
    ----------
    dbdir : str
        Base directory for the database (directory).
    inferemdir : str
        Sub-directory for pickled InferEM instances to evaluate.
    modelprobs : file name
        File name for a pickled dictionary that gives probabilities to use for
        sampling.
    num_sample : int
        The number of machines to sample.
    data : list, None
        Data in a list, or None if prior.
    nprocs : int
        Number of simultaneous processes to run.

    Returns
    -------
    summary_str : str
        A string with computation and timing information.

    """
    from cmpy.inference.bayesianem.util import DictionarySampler

    # initialize list for summary
    summary = []

    # location of inferEM instances
    inferdir = os.path.join(dbdir, inferemdir) 

    # extract information from model probabilities file name
    mp_temp = modelprobs.split('_')
    mp_beta = mp_temp[1]
    mp_penalty = '_'.join([mp_temp[2], mp_temp[3].split('.')[0]])

    # create output directory for samples
    sampledir_name = ''.join(["samples_", 
        inferemdir.split('_')[1],
        "_{}".format(mp_beta),
        "_{}".format(mp_penalty)]
        )
    sampledir = os.path.join(dbdir, sampledir_name)
    if os.path.exists(sampledir):
        raise Exception("Output directory {} exists.".format(sampledir))
    else:
        os.mkdir(sampledir)
        summary.append("*Making directory: {}\n".format(sampledir))

    # add to summary list
    summary.append("* using DB: {}\n".format(dbdir))
    summary.append("  - subdir : {}\n".format(inferdir))
    summary.append("  - samples: {}\n".format(sampledir))

    # open the model probablities dict
    fname = os.path.join(inferdir, modelprobs)
    if os.path.exists(fname):
        summary.append("* Loading {}\n\n".format(fname))
        modelprob_dict = read_probabilities_file(fname)
    else:
        raise Exception("File {} does not exist.".format(fname))

    # intialize dictionary sampler
    machine_name_sampler = DictionarySampler(modelprob_dict)

    # start sampling...
    script_start = datetime.datetime.now()
    script_start_str = script_start.strftime("%H:%M:%S %D")
    summary.append(" -- start time   : {}\n".format(script_start_str))
    
    # create list of sample numbers
    sample_nums = [n for n in range(1, num_sample+1)]
        
    # do the serious computing...
    samples = mp_sample_models(sampledir, sample_nums, machine_name_sampler,
                               data, nprocs)

    # write samples to single file
    f = open(os.path.join(sampledir, 'samples'), 'w')
    # file header
    hdr = ["{}".format('Sample Number'),
           "{}".format('Topology'),
           "{}".format('Type'),
           "{}".format('Cmu'),
           "{}".format('hmu'),
           "{}".format('Edges'),
           "{}".format('States'),
           "{}".format('Start Node'),
           "{}\n".format('Machine Str')]
    f.write("{}".format(','.join(hdr)))
   
    # samples
    for sn in samples:
        f.write("{}".format(','.join(samples[sn])))
    f.close()

    # end processing...
    script_end = datetime.datetime.now()
    script_end_str = script_end.strftime("%H:%M:%S %D")
    summary.append(" -- end time     : {}\n".format(script_end_str))
    time_diff = script_end-script_start
    time_diff_str = deltatime_format(time_diff)
    summary.append(" -- compute time : {}\n\n".format(time_diff_str))

    summary_str = ''.join(summary)

    return summary_str
