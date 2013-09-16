"""util_general.py

Christopher C. Strelioff
chris.stelioff@gmail.com

"""
import os
import argparse
import numpy
import cmpy
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

def create_sample_summary_file(db_dir, sample_dir):
    """Create a file with properties of sample machines from prior or
    posterior.

    Parameters
    ----------
    db_dir : str
        Name of the root database directory.
    sample_dir : str
        Name of the subdirectory containing the sample machine in pickled form.

    """
    # get file 
    sample_files = read_sample_dir(db_dir, sample_dir)

    # open file for writing
    foutname = db_dir + '/' + sample_dir + '.dat'
    fout = open(foutname, 'w')

    # header
    temp = ["Filename",
            "Machine_ID",
            "Recurrent_eM",
            "Excess_Entropy-E",
            "Statistical_Complexity-Cmu",
            "Crypticity-chi",
            "Entropy_Rate-hmu",
            "Gauge_Information-phi",
            "Oracular_Information-zeta",
            "Cryptic_Order",
            "Markov_Order",
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
        
        # estimates from machine, non-L-dependent
        try:
            ee_est = em.excess_entropy()
        except:
            ee_est = numpy.nan
        try:
            cmu_est = em.statistical_complexity()
        except:
            cmu_est = em.stationary_entropy()

        hmu_est = em.entropy_rate()
        
        # number of nodes, orders 
        num_state = len(em.nodes())
        c_ord = em.cryptic_order()
        try:
            m_ord = em.markov_order()
        except:
            # must be a MealyHMM
            m_ord = numpy.nan
        
        # L-dependent values
        cmech=CMechQuantities(em)
        temp = cmech.calculate(['phi', 'zeta', 'chi'], 20)
        # get last (best value) for each quantity
        phi_est = correct_neg_zero(temp[0][-1])  # gauge information
        zeta_est = correct_neg_zero(temp[1][-1]) # oracular information
        chi_est = correct_neg_zero(temp[2][-1])   # crypticity
       
        # get model name and type
        recurr_em = 0
        if type(em) == cmpy.machines.RecurrentEpsilonMachine:
            recurr_em = 1
        em_id = str(em).split(':')[-1].strip()

        temp = ["{fn:s}".format(fn=f),
                "{mn:s}".format(mn=em_id),
                "{id:d}".format(id=recurr_em),
                "{ee:.15e}".format(ee=ee_est),
                "{cmu:.15e}".format(cmu=cmu_est),
                "{chi:.15e}".format(chi=chi_est),
                "{hmu:.15e}".format(hmu=hmu_est),
                "{phi:.15e}".format(phi=phi_est),
                "{zeta:.15e}".format(zeta=zeta_est),
                "{co}".format(co=c_ord),
                "{mo}".format(mo=m_ord),
                "{ns}\n".format(ns=num_state)
            ]
        
        fout.write(' '.join(temp))

    fout.close()

    # change back to startdir
    os.chdir(startdir)

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




