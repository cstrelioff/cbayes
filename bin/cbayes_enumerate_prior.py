#!/usr/bin/env python

"""cbayes_enumerate_prior.py

This script focuses on the prior over models instead of the posterior,
calculating the prior probability for all model topologies in the `machines`
file.

"""
from __future__ import division

import os
import argparse

import cmpy
import cmpy.inference.bayesianem as bayesem

from cbayes import create_machine_prior_file
from cbayes import prior_add_topologies_to_db
from cbayes import check_dir_doesnot_exist

# exception
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
    arg_list.append("-db : Database directory "
            ">> {:s}\n".format(args.database_directory))
    arg_list.append("-nprocs : Number of simultaneous processes to run "
            ">> {:d}\n".format(args.nprocs))
    
    arg_str = ''.join(arg_list)

    return arg_str

def create_parser():
    """Create argparse instance for script.
    
    Returns
    -------
    args : agparse arg instance

    """
    desc_str = (
        """Process all model topologies without data to reflect prior."""
        )
    parser = argparse.ArgumentParser(description=desc_str)
    parser.add_argument('-db', '--database_directory',
            help = 'name of the database directory',
            type = check_dir_doesnot_exist,
            required = True
            )
    parser.add_argument('-nprocs',
            help = 'number of simultaneous processes to run',
            type = int,
            default = 4)
    
    # do the parsing
    args = parser.parse_args()

    return args

def main():
    """Run the specified models, without data, to get InferEM instances
    reflecting the prior.
    
    """
    # parse command line
    args = create_parser()
    
    # get command line args and report settings
    arg_str = report_args(args)
    print arg_str

    # do the serious computing...
    summary_str = create_machine_prior_file(args.database_directory,
                                            args.nprocs)
    # write log
    logfile = os.path.join(args.database_directory, 'summary.log')
    if os.path.exists(logfile):
        f = open(logfile, 'a')
    else:
        f = open(logfile, 'w')
    
    f.write('\n*** start: process prior evidence terms ***\n\n')
    f.write(arg_str)
    f.write('\n')
    f.write(summary_str)
    f.write('\n*** end: process prior evidence terms ***\n')
    f.close()

    print summary_str
    
if __name__ == '__main__':
    # run
    main()

