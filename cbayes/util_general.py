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




