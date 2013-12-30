#!/usr/bin/env python

"""cbayes_enumerate_AddToDB.py

Use the CMPy enumeration method to test all specified uHMM topologies against
the specified datafile.  A pickled dictionary with the evidence (in the
Bayesian sense) for all valid model topologies are written to the specified DB
(really just a directory, to keep things simple) for later use.

"""

from __future__ import division

import os
import argparse

import cmpy
import cmpy.inference.bayesianem as bayesem

from cbayes import read_datafile
from cbayes import check_sr_tuple
from cbayes import add_topologies_to_db

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
    arg_list.append("-a  : Alphabet size >> {:d}\n".format(args.alphabet_size))
    arg_list.append("-n  : Number of states "
            ">> {:s}\n".format(args.number_of_states))
    arg_list.append("--topological_eMs : "
           "topological eMs only? >> {:s}\n".format(str(args.topological_eMs)))
    if args.subsample_range is None:
        arg_list.append("-sr : Subsample range >> "
                "Not provided, *all* data used.\n")
    else:
        pt1, pt2 = args.subsample_range.split(',')
        arg_list.append("-sr : Subsample range >> "
                "{:d}:{:d}\n".format(int(pt1), int(pt2)))
    
    arg_str = ''.join(arg_list)

    return arg_str

def create_parser():
    """Create argparse instance for script.
    
    Returns
    -------
    args : agparse arg instance

    """
    desc_str = (
        """Process candidate model topologies and add to new, or existing,
        DB."""
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
    parser.add_argument('-a', '--alphabet_size',
            help = 'number of letters in alphabet',
            type = int,
            default = 2,
            required = True)
    parser.add_argument('-n', '--number_of_states',
            help = 'a comma sep list of state numbers, eg 1,2,3',
            type = str,
            required = True)
    parser.add_argument('-sr', '--subsample_range',
            help = ("subsample range data[pt1:pt2] passed as -sr pt1,pt2;"
                " *if not provided all data used*"),
            type = check_sr_tuple,
            required = False
            )
    parser.add_argument('--topological_eMs',
            help = 'include only topological eM structures',
            action = 'store_true',
            required = False
            )
    
    # do the parsing
    args = parser.parse_args()

    return args

def main():
    """Run the enumerative inference algorithm with passed data and parameters.
    
    """
    # parse command line
    args = create_parser()

    # process numStates
    ns_str = args.number_of_states.split(',')
    ns_list = [int(n) for n in ns_str]
            
    # read data
    data = read_datafile(args.file)

    # get command line args and report settings
    arg_str = report_args(args)
    print arg_str

    # parse set of models to consider, topological or all?
    if args.topological_eMs:
        em_min = 'min'
    else:
        em_min = 'none'

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
        
    # infer
    # modify directory name to reflect subsampling of data, if any
    (summary_str, inferemdir) = add_topologies_to_db((pt1,pt2), data, 
                                  args.database_directory,
                                  bayesem.LibraryGenerator(args.alphabet_size, 
                                                           ns_list, 
                                                           em_min),
                                  csize=5000)

    # write log
    logfile = os.path.join(inferemdir, 'summary.log')
    if os.path.exists(inferemdir):
        f = open(logfile, 'a')
    else:
        f = open(logfile, 'w')
    
    f.write('\n*** start: Add Models to DB ***\n\n')
    f.write(arg_str)
    f.write('\n')
    f.write(summary_str)
    f.write('\n*** end: Add Models to DB ***\n')
    f.close()

    print summary_str
    
if __name__ == '__main__':
    # run
    main()

