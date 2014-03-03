"""util_general.py

"""
import os
import argparse
import numpy

import cmpy
import cmpy.inference.bayesianem as bayesem
from cmpy.machines.algorithms.mealycmech import CMechQuantities

try:
    import cPickle as pickle
except:
    import pickle

def check_dir_exists(dir):
    """Check if requested output directory exists.

    Paramters
    ---------
    dir : str
        Name of directory to check.

    Exceptions
    ----------
        Raise exception if output directory exists -- avoid overwrite.

    """
    if os.path.exists(dir):
        raise Exception(">> Directory {} exists!".format(dir))

def check_dir_doesnot_exist(dir):
    """Check if requested output directory does not exist.

    Paramters
    ---------
    dir : str
        Name of directory to check.

    Exceptions
    ----------
        Raise exception if output directory does not exist.

    """
    if not os.path.exists(dir):
        raise Exception(">> Directory {} **does not** exist!".format(dir))

    return dir

def check_positive_float(value):
    """Make sure value is positive float."""
    fvalue = float(value)
    
    if fvalue < 0.:
        msg = "{:f} must be a positive float.".format(fvalue)
        raise argparse.ArgumentTypeError(msg)

    return fvalue

def check_positive_int(value):
    """Make sure value is positive int."""
    ivalue = int(value)
    if ivalue <= 0:
        msg = "{:d} must be a positive int.".format(ivalue)
        raise argparse.ArgumentTypeError(msg)

    return ivalue

def check_probability(value):
    """Make sure value is between 0.0 and 1.0."""
    fvalue = float(value)
    if (fvalue < 0.) or (fvalue > 1.):
        msg = "{:f} is not between 0. and 1.".format(fvalue)
        raise argparse.ArgumentTypeError(msg)

    return fvalue

def check_sr_tuple(value):
    "Check sample range tuple."
    try:
        # make sure can split values for sample range
        pt1, pt2 = value.split(',')
        # make sure can cast as int
        ipt1 = int(pt1)
        ipt2 = int(pt2)
    except:
        msg = "argument for -sr must be integer pair; e.g. -sr pt1,pt2"
        raise argparse.ArgumentTypeError(msg)
    
    if ipt1 < 0:
        raise argparse.ArgumentTypeError("{} must be positive".format(ipt1))
    if ipt1 > ipt2:
        raise argparse.ArgumentTypeError("{} must be greater "
                "than {}".format(ipt2, ipt1))

    # return value, unmodified after tests
    return value

def correct_neg_zero(value):
    """Correct values returned as -0.0"""
    if value == -0.0:
        return 0.0
    else:
        return value

def deltatime_format(dt):
    """String format for a datetime.deltatime object.
    
    Parameters
    ----------
    dt : datetime.timedelta instance
        The difference between two datetime.datetime.now() evaluations.
    
    Returns
    dt_str : str
        A formatted string
        
    """
    days = int(dt.days)
    total_seconds = int(dt.total_seconds())
    hours, remainder = divmod(total_seconds,60*60)
    minutes, seconds = divmod(remainder,60)

    dt_str = '{} days {} hrs {} mins {} secs'.format(days,hours,minutes,seconds)

    return dt_str

def create_dir(dir):
    """Create requested output directory.

    Parameters
    ----------
    dir : str
        Name of directory to create in current working directory.

    Returns
    -------
    dir : str
        Path to created output directory.

    Exceptions
    ----------
        Raise exception if output directory exists -- avoid overwrite.

    """
    cwd = os.getcwd()
    outdir = os.path.join(cwd, dir)

    if os.path.exists(outdir):
        raise Exception(">> Directory {} exists!".format(outdir))
    else:
        print ("Making output directory: {}\n".format(outdir))
        os.mkdir(outdir)

    return outdir

def read_datafile(filename):
    """Read data file.

    Parameters
    ----------
    filename : str
        Data file name -- assumed to contain comma-separated series of 0s, 1s,
        etc.  Lines starting with `#` are ignored.

    Returns
    -------
    data : list
        A list containing the data series.

    """
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                pass
            else:
                line = line.strip()
                data += line.split(',')

    return data

def read_evidence_file(inferemdir):
    """Read evidence file in passed `inferdir` directory.

    Parameters
    ----------
    inferemdir : str
        Full path to directory containing the `evidence` file.

    Returns
    -------
    evidence_dict : dict
        A dictionary of evidence values from specified file.

    """
    # open machines file and read information
    efile = os.path.join(inferemdir, 'log_evidence')
    f = open(efile, 'r')
    evidence_raw = f.readlines()
    f.close()

    # process and return
    headerline = True
    evidence = {}
    for line in evidence_raw:
        if line.startswith('#'):
            pass
        else:
            if headerline:
                # hit header row, next row real data
                headerline = False
            else:
                line = line.strip()
                # em_name, log_evidence, nodes, edges
                em_name, log_evidence, nodes, edges = line.split(',')
                evidence[em_name] = {'log_evidence': float(log_evidence),
                                     'nodes': int(nodes),
                                     'edges': int(edges)}

    return evidence

def read_machines_file(db):
    """Read machines file in passed `db` directory.

    Parameters
    ----------
    db : os.path instance
        Full path to directory containing the `machines` file.
        
    Returns
    -------
    machines : list
        A list of lists, containg information on all candidate machine
        topologies.

    """
    # open machines file and read information
    mfile = os.path.join(db, 'machines')
    f = open(mfile, 'r')
    machines_raw = f.readlines()
    f.close()

    # process and return
    headerline = True
    machines = []
    for line in machines_raw:
        if line.startswith('#'):
            pass
        else:
            if headerline:
                # hit header row, next row real data
                headerline = False
            else:
                line = line.strip()
                # em_name, em_type, em_num_states, em_num_edges, em_str
                machines.append(line.split(','))
    
    return machines

def read_probabilities_file(filename):
    """Read a file containing the machine/model probabilities.

    Parameters
    ----------
    filename : str
        Name of the probabilities file -- complete relative path

    Returns
    -------
    probabilities : dict
        A dictionary of model probabilties

    """
    # open file
    f = open(filename, 'r')
    probs_raw = f.readlines()
    f.close()

    # process and return
    headerline = True
    probabilities = {}
    for line in probs_raw:
        if line.startswith('#'):
            pass
        else:
            if headerline:
                # hit header row, next row real data
                headerline = False
            else:
                line = line.strip()
                em_name, prob = line.split(',')

                # save machine probs to dict
                probabilities[em_name] = float(prob)

    return probabilities

def read_sample_dir(db_dir, sample_dir):
    """Read sample files from given directory.
    
    Parameters
    ----------
    db_dir : str
        Name of root database directory.
    sample_dir : str
        Name of sub-directory that contains the sampled machines.
    
    """
    # save startdir
    startdir = os.getcwd()

    # specify set of samples using dir name
    os.chdir(os.path.join(db_dir, sample_dir))
    cwd = os.getcwd()

    # read files in cwd
    files = os.listdir('.')

    # make sure all files are sample files and sort before return
    def check_sample_file(filename):
        if filename.startswith('sample') and filename.endswith('pickle'):
            return True
        else:
            return Fale

    sample_files = [f for f in sorted(files) if check_sample_file(f)]

    # change directory back to start dir
    os.chdir(startdir)

    return sample_files

def write_data(directory, file, process, length, eM):
    """Write the data file.
    
    Parameters
    ---------
    directory: str
        Directory to write file -- '.' is a good option.
    file : str
        Assumed to be name of data file, for example Even.dat.
    process : str
        Name of the process.
    length : int
        Length of the desired data series.
    eM : RecurrentEpsilonMachine, MealyHMM
        Machine to use for data generation.
    
    """
    
    cwd = os.getcwd()
    out_dir = cwd
    if directory != '.':
        out_dir = os.path.join(cwd, directory)
        if os.path.exists(out_dir):
            os.chdir(out_dir)
        else:
            raise ProcessException("Output directory does not exist.")

    f = open(file, 'w')

    # header
    f.write("# SETTINGS:\n")
    f.write("# -f    : Output data file >> {:s}\n".format(file))
    f.write("# -l    : Length of data series >> {:d}\n".format(length))
    f.write("# -p    : Process >> {:s}\n".format(process))
    f.write("# -d    : Ouput diretory >> {:s}\n".format(directory))
    f.write('#\n')

    # eM details
    f.write('# eM \n')
    f.write("# name: {:s}\n".format(str(eM)))

    trans_str = eM.to_string()
    trans_list = trans_str.split(';')
    for element in trans_list:
        element = element.strip()
        f.write('# ' + element + '\n')

    # now, data
    eM.randomize_current_node()
    sN = eM.get_current_node()
    f.write("# start node: {:s}\n".format(str(sN)))
    f.write('#\n')
    
    # create the data
    data = eM.symbols(length)
    data_str = ','.join(data)

    # find the evidence for generated data
    evidence_available = True
    try:
        posterior = bayesem.InferEM(eM, data)
        log_evidence = posterior.log_evidence()
        for beta in [0.0, 2.0, 4.0]:
            curr_evi = log_evidence - beta*len(eM.nodes())
            tmp= ("# log P(D|Mi)-beta*|Mi| "
                "(beta={b:f}) = {ce:.4f}\n"
                )
            f.write(tmp.format(b=beta, ce=curr_evi))
        f.write("#\n")
    except:
        evidence_available = False

    # compare with IID
    iid = cmpy.machines.IID(len(eM.alphabet()))
    iid_posterior = bayesem.InferEM(iid, data)
    iid_evidence = iid_posterior.log_evidence()
    f.write("# log P(D|IID) (beta=0.0) = {:.4f}\n".format(iid_evidence))
    f.write("#\n")

    # write data and close file
    f.write(data_str)
    f.close()

    # include an inference log file
    if evidence_available:
        flog = open(file.split('.')[0] + '.log', 'w')
        flog.write("# Inference information for {:s}\n".format(file))
        flog.write("Evi_True_beta_0 Evi_True_beta_2 ")
        flog.write("Evi_True_beta_4 Evi_IID_beta_0\n")
        flog.write("{:.2f} ".format(log_evidence))
        flog.write("{:.2f} ".format(log_evidence - 2.*len(eM.nodes())))
        flog.write("{:.2f} ".format(log_evidence - 4.*len(eM.nodes())))
        flog.write("{:.2f}\n".format(iid_evidence))
        flog.close()
    else:
        flog = open(file.split('.')[0] + '.log', 'w')
        flog.write("# Inference information for {:s}\n".format(file))
        flog.write("Evi_IID_beta_0\n")
        flog.write("{:.2f}\n".format(iid_evidence))
        flog.close()

    # change back to entry directory
    os.chdir(cwd)

def write_em_pickle(file, eM):
    """Write a pickle file of the generated process for later use.
    
    Parameters
    ---------
    file : str
        Assumed to be name of data file, for example Even.dat
    eM : RecurrentEpsilonMachine, MealyHMM
        Machine to pickle.

    """
    try:
        filename = file.split('.')[0]
    except:
        filename = file
    
    f = open(filename + '.pickle', 'w')
    pickle.dump(eM, f)
    f.close()

def write_evidence_file(evidence, filename):
    """Write a file containing the machine/model evidence terms.

    Parameters
    ----------
    evidence : dict
        A dictionary containg the evidence data
    filename : str
        Name for output file

    Output
    ------
    log_evidence : file
        Write log_evidence file

    """
    # open file
    f = open(filename, 'w')
    
    # write header
    f.write('{},{},{},{}\n'.format('em_name', 'log_evidence',
                                   'nodes', 'edges'))

    for em_name in sorted(evidence.iterkeys()):
        em_info = evidence[em_name]
        f.write('{},{},{},{}\n'.format(em_name,
                                     em_info['log_evidence'],
                                     em_info['nodes'],
                                     em_info['edges']))
    f.close()

def write_probabilities_file(probabilities, filename):
    """Write a file containing the machine/model probabilities.

    Parameters
    ----------
    probabilities : dict
        A dictionary containg the probabilities for models
    filename : str
        Name for output file

    Output
    ------
    mpfile : file
        Write model probability file

    """
    # open file
    f = open(filename, 'w')

    # write header
    f.write('{},{}\n'.format('em_name', 'probability'))

    for em_name in sorted(probabilities.iterkeys()):
        f.write('{},{}\n'.format(em_name, probabilities[em_name]))

    f.close()
