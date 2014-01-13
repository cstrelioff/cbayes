#!/usr/bin/env python

"""cbayes_enumerate_ProcessSamples.py

This script processes the samples generated using `cbayes_enumerate_Sample.py`.
A single *.dat file is created with basic information about each sampled
machine including machine id, number of states, number of edges, Cmu, hmu etc.
--- see `cbayes.util_infer_db.py` for details.

"""
from __future__ import division

import os
import argparse

import cmpy
import cmpy.inference.bayesianem as bayesem

from cbayes import check_positive_float
from cbayes import create_sample_summary_file

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
    arg_list.append("-db : Database root directory "
            ">> {:s}\n".format(args.database_directory))
    arg_list.append("-sdir : Sample sub-directory "
            ">> {:s}\n".format(args.sample_directory))
    
    arg_str = ''.join(arg_list)

    return arg_str

def create_parser():
    """Create argparse instance for script.
    
    Returns
    -------
    args : agparse arg instance

    """
    desc_str = (
        """Create a summary file from a directory of sampled machines using
        `cbayes_enumerate_Sample.py`."""
        )
    parser = argparse.ArgumentParser(description=desc_str)
    parser.add_argument('-db', '--database_directory',
            help = 'name of the database directory',
            type = str,
            required = True
            )
    parser.add_argument('-sdir', '--sample_directory',
            help = 'name of subdirectory with desired sample (machine) pickles',
            type = str,
            required = True
            )
    
    # do the parsing
    args = parser.parse_args()

    return args

def main():
    """Produce a summary data file for a directory of pickled
    `RecurentEpsilonMachines` or `MeallyHMMs` that represent samples from the
    prior or posterior distirbutions. Create these using available scripts:

    * `cbayes_enumerate_AddToDB.py` or `cbayes_enumerate_PriorAddToDB.py`
    * `cbayes_enumerate_CalcProb.py`
    * `cbayes_enumerate_Sample.py`
    
    """
    # parse command line
    args = create_parser()

    # get command line args and report settings
    arg_str = report_args(args)
    print arg_str
    
    # create sample summary file
    summary_str = create_sample_summary_file(args.database_directory,
                                             args.sample_directory)


    ## print compute time summary
    print summary_str
    
if __name__ == '__main__':
    main()

