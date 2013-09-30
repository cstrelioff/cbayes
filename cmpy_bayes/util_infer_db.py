"""util_infer_db.py

Christopher C. Strelioff
chris.strelioff@gmail.com

"""
from __future__ import absolute_import
from __future__ import division

import numpy
from multiprocessing import Pool, cpu_count
import itertools
import datetime
import os
from math import exp
try:
    import cPickle as pickle
except:
    import pickle

import cmpy
from cmpy.math import log, logaddexp
import cmpy.inference.bayesianem as cmbayes
import cmpy.orderlygen.pyicdfa as pyidcdfa

from .util_general import read_sample_dir

def create_sample_summary_file(db_dir, sample_dir):
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
    script_start = datetime.datetime.now()
    summary = []
    summary.append(" -- start time   : {}\n".format(script_start))

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
    fout.write(' '.join(temp))

    # change dir
    startdir = os.getcwd()
    os.chdir(os.path.join(db_dir, sample_dir))

    ## cycle through files
    for f in sample_files:
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
        
        fout.write(' '.join(temp))

    fout.close()

    # change back to startdir
    os.chdir(startdir)
    
    script_end = datetime.datetime.now()
    summary.append(" -- end time     : {}\n".format(script_end))
    time_diff = str(script_end-script_start)
    summary.append(" -- compute time : {}\n\n".format(time_diff))

    summary_str = ''.join(summary)

    return summary_str

def infer_map(machine, data):
    """A map function for mutiprocessing to consider topologies in parallel.
    
    Return:
    
    (mname, minfo): tuple
        Tuple with machine name `mname` and dictionary `minfo` containing
        `log_evidence` for data series provided as well as the number of
        `nodes` and `edges` for the machine.
    
    """
    # pass to InferEM
    infer_temp = cmbayes.InferEM(machine,data)

    # extract machine name and other information
    model_info = (str(machine),{'log_evidence': infer_temp.log_evidence(),
                                'nodes': len(machine.nodes()),
                                'edges': len(machine.edges())})
    
    return model_info

def infer_map_star(a_b):
    """Convert `infer_map([1,2])` to `infer_map(1,2)` call."""
    return infer_map(*a_b)

def model_probs_factor(inferemfile, args):
    """Calculate model probability factor for model given beta 
    and penalty type.
    """
    beta = args[0]
    penalty = args[1]
    
    f = open(inferemfile, 'rb')
    temp = pickle.load(f)
    f.close()
    
    # get evidence
    curr_evidence = temp.log_evidence()

    # modify with beta and penalty
    penalty_term = 0.
    if penalty == 'num_states':
        # get number of nodes
        penalty_term = len(temp.machine.nodes())
    elif penalty == 'num_edges':
        # get number of edges
        penalty_term = len(temp.machine.edges())

    # for each machine topology, add log evidence
    evidence = curr_evidence - beta*(penalty_term)
    
    # construct short filename
    shortfilename = inferemfile.split('/')[-1]
    
    # return filename, evidence
    return (shortfilename, evidence)
    
def model_probs_factor_star(a_b):
    """Convert `model_probs_factor([1,2])` to `model_probs_factor(1,2)`
    call.
    
    """
    
    return model_probs_factor(*a_b)

def sample_map(sample_num, args):
    """A map function for mutiprocessing to consider sampling of machines.
    
    Parameters
    ----------
    sample_num : int
        The sample number -- used for naming output files uniquely.
    args :
        Used to pass sampler and outdir.

    """
    data = args[0]
    sampler = args[1]
    inputdir = args[2]
    outdir = args[3]

    # create file name
    fname = "sample{:05d}.pickle".format(sample_num)

    # sample a modelname
    sample_mname = sampler()
    
    # create filename for sampled model name
    sample_file = ''.join(['inferEM_', sample_mname, '.pickle'])
    inferem_file = os.path.join(inputdir, sample_file)

    # if exists, load.  otherwise create
    if os.path.exists(inferem_file):
        # pickled inferEM exists
        f = open(inferem_file, 'r')
        inferem_instance = pickle.load(f)
        f.close()
    else:
        # pickled inferEM does not exist, create
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
        inferem_instance = cmbayes.InferEM(machine, data)

    # sample machine (also returns start node, not needed)
    _, em_sample = inferem_instance.generate_sample()
    
    # write to file
    f = open(os.path.join(outdir,fname), 'w')
    
    pickle.dump(em_sample, f)
    f.close()

    return sample_num

def sample_map_star(a_b):
    """Convert `sample_map_star([1,2,3])` to `sample_map(1,2,3)` call."""
    
    return sample_map(*a_b)

def add_topologies_to_db(range, data=None, dbdir=None, iter_topologies=None,
        csize=1000):
    """A fuction to add new topologies to a database consisting of a pickled
    dictionary with model id and evidence value for provided data.  If the
    passed directory does not exist, it will be created.

    Parameters
    ----------
    range : (int, int) tuple
        Subsample range for the passed data.  That is, the data analyzed is
        already subsampled to the range d[pt1:pt2] if (pt1,pt2) is passed.
    data : list
        A list containg the data to be used for inference.
    dbdir : str
        The database (directory) name.
    iter_topologies : iterator of MealyHMM or RecurrentEpsilonMachines
        An iterator with candidate topologies.
    csize : int
        Chunk size for the multi-threading.

    Returns
    -------
    (summary_str, inferdir) : (str, str)
        A summary of the inference done for print to stdout or file and the
        path to the InferEM subdirectory.

    """
    cwd = os.getcwd()
    summary = []

    # check data
    if (data == None) or (dbdir == None) or (iter_topologies == None):
        raise Exception("""
        Must pass data, database directory and iterator of
        topologies to be considered.
        """)
    
    # check if database dir exisits
    inferdir = os.path.join(dbdir,"inferEM_{:d}-{:d}".format(range[0], 
                                                             range[1]))
    if os.path.exists(dbdir):
        summary.append("* DB {:s} exists\n".format(dbdir))
        # test for inferEM directory
        if os.path.exists(inferdir):
            # okay, pass
            pass 
        else:
            summary.append("* making subdirectory: {}\n".format(inferdir))
            os.mkdir(inferdir)
    else:
        # if not, create
        summary.append("* making DB: {}\n".format(dbdir))
        os.mkdir(dbdir)
        summary.append("* making subdirectory: {}\n".format(inferdir))
        os.mkdir(inferdir)
    
    # change to inferEM database directory
    os.chdir(inferdir)
    summary.append("* using DB: {} subdir: {}\n".format(dbdir, inferdir))

    # find number of processes possible on machine
    numCPUs = cpu_count()
    summary.append("* Starting pool with {} processes\n\n".format(numCPUs))
    # start a pool with this number of processes
    pool = Pool(processes=numCPUs)
    
    # start inference...
    script_start = datetime.datetime.now()
    summary.append(" -- start time   : {}\n".format(script_start))

    # run with imap to get iterator execution
    evidence_dict = {}
    num_tried = 0
    num_valid = 0

    #
    # use imap with some manipulation to take multiple arguments
    #  http://stackoverflow.com/questions/5442910/
    #         python-multiprocessing-pool-map-for-multiple-arguments
    #
    for mname,minfo in pool.imap(infer_map_star, itertools.izip(iter_topologies,
                                                       itertools.repeat(data)),
                                                       chunksize=csize):
        # count number tried
        num_tried += 1
        
        # m_info format: {'id': str, 'log_evidence': float}
        if minfo['log_evidence'] == -numpy.inf:
            # no valid for this data and machine
            pass
        else:
            num_valid += 1
            evidence_dict[mname] = minfo

    # check for existing dicionary
    num_existed = 0
    if os.path.isfile('evidence_dictionary.pickle'):
        summary.append(" ** evidence_dictionary.pickle already exists!\n")
        summary.append("   ** only adding new models.\n")
        summary.append("   ** existing enteries will NOT be modified.\n")

        # load existing file
        f = open('evidence_dictionary.pickle', 'rb')
        evidence_previous = pickle.load(f)
        f.close()

        # parse new results
        for mid, mid_evid in evidence_dict.items():
            if evidence_previous.has_key(mid):
                # already exists
                num_existed += 1
            else:
                # new entry, add
                evidence_previous[mid] = mid_evid

        # pickle the evidence dictionary, new and old results
        f = open('evidence_dictionary.pickle', 'wb')
        pickle.dump(evidence_previous,f)
        f.close()
    else:
        # pickle the evidence dictionary, only new results
        f = open('evidence_dictionary.pickle', 'wb')
        pickle.dump(evidence_dict,f)
        f.close()
    
    script_end = datetime.datetime.now()
    summary.append(" -- end time     : {}\n".format(script_end))
    time_diff = str(script_end-script_start)
    summary.append(" -- compute time : {}\n\n".format(time_diff))
    summary.append(" --  {} valid of {} attempted "
            "({} already existed)\n".format(num_valid, num_tried, num_existed))
    
    summary_str = ''.join(summary)

    # return to cwd
    os.chdir(cwd)

    return (summary_str, inferdir)

def calc_probs_beta_db(dbdir, inferemdir, beta, penalty):
    """A function to calculate the probabilities for a pickled dictionary of
    model evidence values.

    Parameters
    ----------

    dbdir : str
        Base directory for the database (directory).
    inferemdir : str
        Sub-directory for pickled InferEM instances to evaluate.
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

    # check if database dir exisits
    if os.path.exists(dbdir):
        # test for inferEM directory
        if os.path.exists(inferdir):
            # okay, pass
            pass 
        else:
            exit('\n\n ** InferEM subdirectory not found in %s !\n' % dbdir)
    else:
        # if not, error
        exit('\n\n ** Database directory %s not found!\n' % dbdir)

    # check for existing output file before doing work
    ftemp = "modelprobs_beta-{:f}_penalty-{:s}.pickle".format(beta, penalty)
    fname = os.path.join(inferdir, ftemp)
    if os.path.exists(fname):
        raise Exception("\nOutput file {} exists!\n".format(fname))

    # add to summary list
    summary.append("* using DB: {}\n".format(dbdir))
    summary.append("  - subdir: {}\n\n".format(inferdir))
    
    # start inference...
    script_start = datetime.datetime.now()
    summary.append(" -- start time   : {}\n".format(script_start))

    # open evidence dictionary from disk
    fevi = open(os.path.join(inferdir, "evidence_dictionary.pickle"), 'rb')
    evidence_dict = pickle.load(fevi)
    fevi.close()
    
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
        model_probabilities[em] = exp(temp_log_prob)
    
    # save pickled instance of dictionary
    f = open(fname, 'w')
    pickle.dump(model_probabilities,f)
    f.close()

    script_end = datetime.datetime.now()
    summary.append(" -- end time     : {}\n".format(script_end))
    time_diff = script_end - script_start
    summary.append(" -- compute time : {}\n".format(str(time_diff)))

    summary_str = ''.join(summary)

    return (summary_str, model_probabilities)


def sample_db(data, dbdir, inferemdir, modelprobs, num_sample, prior=False):
    """A function to generate samples from the prior or posterior over model
    topologies.

    Parameters
    ----------
    data : list
        The data used for inference.
    dbdir : str
        Base directory for the database (directory).
    inferemdir : str
        Sub-directory for pickled InferEM instances to evaluate.
    modelprobs : file name
        File name for a pickled dictionary that gives probabilities to use for
        sampling.
    num_sample : int
        The number of machines to sample.
    prior : bool
        Indicate this is sampling from prior -- ignore data.
    
    Returns
    -------

    """
    from cmpy.inference.bayesianem.util import DictionarySampler

    # initialize list for summary
    summary = []

    # location of inferEM instances
    inferdir = os.path.join(dbdir, inferemdir) 

    # check if database dir exisits
    if os.path.exists(dbdir):
        # test for inferEM directory
        if os.path.exists(inferdir):
            # okay, pass
            pass 
        else:
            exit('\n\n ** InferEM subdirectory not found in %s !\n' % dbdir)
    else:
        # if not, error
        exit('\n\n ** Database directory %s not found!\n' % dbdir)

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
        #summary.append("*Adding to directory: {}\n".format(sampledir))
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
        f = open(fname, 'r')
        modelprob_dict = pickle.load(f)
        f.close()
    else:
        raise Exception("File {} does not exist.".format(fname))

    # intialize dictionary sampler
    file_sampler = DictionarySampler(modelprob_dict)

    # find number of processes possible on machine
    numCPUs = cpu_count()
    summary.append("* Starting pool with {} processes\n".format(numCPUs))

    # start a pool with this number of processes
    pool = Pool(processes=numCPUs)
    
    # start inference...
    script_start = datetime.datetime.now()
    summary.append(" -- start time   : {}\n".format(script_start))
    
    # define iterator fuction for sample number
    def iter_samplenum(sample_num):
        for n in range(sample_num):
            yield n
    
    # create instance of iter
    iter_snum = iter_samplenum(num_sample)
    
    #
    # use imap with some manipulation to take multiple arguments
    #  http://stackoverflow.com/questions/5442910/
    #         python-multiprocessing-pool-map-for-multiple-arguments
    #
    # limit chunksize to 1000
    csize = min(1000, num_sample/numCPUs)
    if prior:
        # this is prior, pass data=None
        args = [None, file_sampler, inferdir, sampledir]
    else:
        # posterior, pass real data
        args = [data, file_sampler, inferdir, sampledir]

    for snum in pool.imap(sample_map_star, 
                          itertools.izip(iter_snum,
                                         itertools.repeat(args)), 
                          chunksize=csize):
        # do nothing
        pass
    
    pool.close()
    pool.join()

    script_end = datetime.datetime.now()
    summary.append(" -- end time     : {}\n".format(script_end))
    time_diff = script_end - script_start
    summary.append(" -- compute time : {}\n".format(str(time_diff)))

    summary_str = ''.join(summary)

    return summary_str

def prior_add_topologies_to_db(dbdir=None, iter_topologies=None, csize=1000):
    """A fuction to add new topologies to a database consisting of a pickled
    dictionary with model id and evidence value for provided data.  If the
    passed directory does not exist, it will be created.

    Parameters
    ----------
    dbdir : str
        The database (directory) name.
    iter_topologies : iterator of MealyHMM or RecurrentEpsilonMachines
        An iterator with candidate topologies.
    csize : int
        Chunk size for the multi-threading.

    Returns
    -------
    (summary_str, inferdir) : (str, str)
        A summary of the inference done for print to stdout or file and the
        path to the InferEM subdirectory.

    """
    cwd = os.getcwd()
    summary = []

    # check data
    if (dbdir == None) or (iter_topologies == None):
        raise Exception("""
        Must pass database directory and iterator of topologies to be
        considered.
        """)
    
    # check if database dir exisits
    inferdir = os.path.join(dbdir,"inferEM_{:d}-{:d}".format(0, 0))
    if os.path.exists(dbdir):
        summary.append("* DB {:s} exists\n".format(dbdir))
        # test for inferEM directory
        if os.path.exists(inferdir):
            # okay, pass
            pass 
        else:
            summary.append("* making subdirectory: {}\n".format(inferdir))
            os.mkdir(inferdir)
    else:
        # if not, create
        summary.append("* making DB: {}\n".format(dbdir))
        os.mkdir(dbdir)
        summary.append("* making subdirectory: {}\n".format(inferdir))
        os.mkdir(inferdir)
    
    # change to inferEM database directory
    os.chdir(inferdir)
    summary.append("* using DB: {} subdir: {}\n".format(dbdir, inferdir))

    # find number of processes possible on machine
    numCPUs = cpu_count()
    summary.append("* Starting pool with {} processes\n\n".format(numCPUs))
    # start a pool with this number of processes
    pool = Pool(processes=numCPUs)
    
    # start inference...
    script_start = datetime.datetime.now()
    summary.append(" -- start time   : {}\n".format(script_start))

    # run with imap to get iterator execution
    evidence_dict = {}
    num_tried = 0
    num_valid = 0

    #
    # use imap with some manipulation to take multiple arguments
    #  http://stackoverflow.com/questions/5442910/
    #         python-multiprocessing-pool-map-for-multiple-arguments
    #
    ## pass None as data to get prior
    for mname,minfo in pool.imap(infer_map_star, itertools.izip(iter_topologies,
                                                       itertools.repeat(None)),
                                                       chunksize=csize):
        # count number tried
        num_tried += 1
        
        # m_info format: {'id': str, 'log_evidence': float}
        if minfo['log_evidence'] == -numpy.inf:
            # not valide -- should not visit this for prior...
            pass
        else:
            num_valid += 1
            evidence_dict[mname] = minfo

    # check for existing dicionary
    num_existed = 0
    if os.path.isfile('evidence_dictionary.pickle'):
        summary.append(" ** evidence_dictionary.pickle already exists!\n")
        summary.append("   ** only adding new models.\n")
        summary.append("   ** existing enteries will NOT be modified.\n")

        # load existing file
        f = open('evidence_dictionary.pickle', 'rb')
        evidence_previous = pickle.load(f)
        f.close()

        # parse new results
        for mid, mid_evid in evidence_dict.items():
            if evidence_previous.has_key(mid):
                # already exists
                num_existed += 1
            else:
                # new entry, add
                evidence_previous[mid] = mid_evid

        # pickle the evidence dictionary, new and old results
        f = open('evidence_dictionary.pickle', 'wb')
        pickle.dump(evidence_previous,f)
        f.close()
    else:
        # pickle the evidence dictionary, only new results
        f = open('evidence_dictionary.pickle', 'wb')
        pickle.dump(evidence_dict,f)
        f.close()
    
    script_end = datetime.datetime.now()
    summary.append(" -- end time     : {}\n".format(script_end))
    time_diff = str(script_end-script_start)
    summary.append(" -- compute time : {}\n\n".format(time_diff))
    summary.append(" --  {} valid of {} attempted "
            "({} already existed)\n".format(num_valid, num_tried, num_existed))
    
    summary_str = ''.join(summary)

    # return to cwd
    os.chdir(cwd)

    return (summary_str, inferdir)

