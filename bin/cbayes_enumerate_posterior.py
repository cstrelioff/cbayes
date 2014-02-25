#!/usr/bin/env python

"""cbayes_enumerate_posterior.py

This script focuses on the posterior over models, calculating the posterior log
evidence for all model topologies in the `machines` file using the specified
data file and data range.

"""

from __future__ import division

import os
import argparse

import cmpy
import cmpy.inference.bayesianem as bayesem

from cbayes import check_sr_tuple
from cbayes import create_machine_posterior_file
from cbayes import read_datafile

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
    arg_list.append("-f  : Data file >> {:s}\n".format(args.file))
    arg_list.append("-db : Database directory "
            ">> {:s}\n".format(args.database_directory))
    if args.subsample_range is None:
        arg_list.append("-sr : Subsample range >> "
                "Not provided, *all* data used.\n")
    else:
        pt1, pt2 = args.subsample_range.split(',')
        arg_list.append("-sr : Subsample range >> "
                "{:d}:{:d}\n".format(int(pt1), int(pt2)))
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
        """Use machines file in specified db to do inference on data file
        specified by --file."""
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
    parser.add_argument('-sr', '--subsample_range',
            help = ("subsample range data[pt1:pt2] passed as -sr pt1,pt2;"
                " *if not provided all data used*"),
            type = check_sr_tuple,
            required = False
            )
    parser.add_argument('-nprocs',
            help = 'number of simultaneous processes to run',
            type = int,
            default = 4)
    
    
    # do the parsing
    args = parser.parse_args()

    return args

def main():
    """Run the enumerative inference algorithm with passed data and specified
    machines file.
    
    """
    # parse command line
    args = create_parser()

    # read data
    data = read_datafile(args.file)

    # get command line args and report settings
    arg_str = report_args(args)
    print arg_str

    # process subsample range, if any
    if args.subsample_range is None:
        # use all data
        pt1 = 0
        pt2 = len(data)
    else:
        # use part of data
        pt1, pt2 = args.subsample_range.split(',')
        pt1 = int(pt1)
        pt2 = int(pt2)
        data = data[pt1:pt2]
    
    # do the serious computing...
    summary_str = create_machine_posterior_file(args.database_directory,
                                                (pt1, pt2),
                                                data,
                                                args.nprocs)
    
    # write log
    logfile = os.path.join(args.database_directory, 'summary.log')
    if os.path.exists(logfile):
        f = open(logfile, 'a')
    else:
        f = open(logfile, 'w')
    
    f.write("""\n*** start: process posterior evidence terms,"""
            """ range {}, {} ***\n\n""".format(pt1, pt2))
    f.write(arg_str)
    f.write('\n')
    f.write(summary_str)
    f.write("""\n*** end: process posterior evidence terms,"""
            """ range {}, {} ***\n""".format(pt1, pt2))
    f.close()

    print summary_str

if __name__ == '__main__':
    main()

