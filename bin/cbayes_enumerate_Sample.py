#!/usr/bin/env python

"""cbayes_enumerate_Sample.py

This script generates samples from the prior, or posterior, over model
topologies using the output from other scripts.  Create pickled InferEM
instances using

* `cbayes_enumerate_AddToDB.py`
or
* `cbayes_enumerate_PriorAddToDB.py`

A dictionary of model probabilities can then be created using output from the
above scripts and running

* `cbayes_enumerate_CalcProbs.py`

This script uses these outputs to generate sample epsilon-machines or uHMMs and
writes pickled instances of the machines to a directory in the DB.

"""
from __future__ import division

import os
import argparse

import cmpy
import cmpy.inference.bayesianem as bayesem

from cmpy_bayes import read_datafile
from cmpy_bayes import check_positive_float
from cmpy_bayes import sample_db

class InferDBException(Exception):
    pass

def report_args(args):
    """Report the requested settings.
    
    Returns
    -------
    arg_str : str
        A summary of the script settings.

    """
    arg_list =[]

    arg_list.append("SETTINGS:\n")
    arg_list.append("-f  : Data file >> {:s}\n".format(args.file))
    arg_list.append("-db : Database root directory "
            ">> {:s}\n".format(args.database_directory))
    arg_list.append("-idir : InferEM  sub-directory "
            ">> {:s}\n".format(args.inferem_directory))
    arg_list.append("-mp : Model probabilities "
            ">> {:s}\n".format(args.model_probabilities))
    arg_list.append("-ns : Number of machines to sample "
            ">> {:d}\n".format(args.number_samples))
    
    arg_str = ''.join(arg_list)

    return arg_str

def create_parser():
    """Create argparse instance for script.
    
    Returns
    -------
    args : agparse arg instance

    """
    desc_str = (
        """Sample machines using a directory of pickled InferEM instances and 
        dictionary of model probabiltiies."""
        )
    parser = argparse.ArgumentParser(description=desc_str)
    parser.add_argument('-f','--file', 
            help = 'input (data) file name', 
            type = str,
            required = True)
    parser.add_argument('-db', '--database_directory',
            help = 'name of the database directory',
            type = str,
            required = True
            )
    parser.add_argument('-idir', '--inferem_directory',
            help = 'name of subdirectory with desired InferEM pickles',
            type = str,
            required = True
            )
    parser.add_argument('-mp', '--model_probabilities',
            help = 'name of file with pickled dictionary of probabilities',
            type = str,
            required = True
            )
    parser.add_argument('-ns', '--number_samples',
            help = 'number of machines to sample',
            type = int,
            required = True
            )
    
    # do the parsing
    args = parser.parse_args()

    return args

def main():
    """Sample machines using a directory of pickled InferEM files and a pickled
    dictionary of model probabilities.  Create these using available scripts:

    * cbayes_enumerate_AddToDB.py
    * cbayes_enumerate_CalcProb.py
    
    """
    # parse command line
    args = create_parser()

    # get command line args and report settings
    arg_str = report_args(args)
    print arg_str

    # read data
    data = read_datafile(args.file)
    
    # call sample_db method
    summary_str = sample_db(data, args.database_directory,
            args.inferem_directory, args.model_probabilities,
            args.number_samples)
    
    print summary_str

    # write log
    inferemdir = os.path.join(args.database_directory, args.inferem_directory)
    logfile = os.path.join(inferemdir, 'summary.log')
    if os.path.exists(inferemdir):
        f = open(logfile, 'a')
    else:
        f = open(logfile, 'w')
     
    f.write('\n*** start: Generate Sample Machines***\n\n')
    f.write(arg_str)
    f.write('\n')
    f.write(summary_str)
    f.write('\n*** end: Generate Sample Machines***\n')
    f.close()
    
if __name__ == '__main__':
    main()

