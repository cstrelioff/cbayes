#!/usr/bin/env python

"""cbayes_enumerate_PriorAddToDB.py

This script is similar to cbayes_enumerate_AddToDB.py except that the focus is
on the prior over models instead of the posterior.  The objective is to
understand what the setting for hyper-parameters of priors at various levels
say about stated a priori assumptions.  The resulting directory can be used for
sampling from the prior over model topologies.

"""
from __future__ import division

import os
import argparse

import cmpy
import cmpy.inference.bayesianem as bayesem

from cmpy_bayes import prior_add_topologies_to_db

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
    arg_list.append("-a  : Alphabet size >> {:d}\n".format(args.alphabet_size))
    arg_list.append("-n  : Number of states "
            ">> {:s}\n".format(args.number_of_states))
    arg_list.append("--topological_eMs : "
           "topological eMs only? >> {:s}\n".format(str(args.topological_eMs)))
    
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
    parser.add_argument('--topological_eMs',
            help = 'include only topological eM structures',
            action = 'store_true',
            required = False
            )
    
    # do the parsing
    args = parser.parse_args()

    return args

def main():
    """Run the specified models, without data, to get InferEM instances
    reflecting the prior.
    
    """
    # parse command line
    args = create_parser()

    # process numStates
    ns_str = args.number_of_states.split(',')
    ns_list = [int(n) for n in ns_str]
    
    # get command line args and report settings
    arg_str = report_args(args)
    print arg_str

    # parse set of models to consider, topological or all?
    if args.topological_eMs:
        em_min = 'min'
    else:
        em_min = 'none'

    # infer
    (summary_str, inferemdir) = prior_add_topologies_to_db(\
                                  args.database_directory,
                                  bayesem.LibraryGenerator(args.alphabet_size, 
                                                           ns_list, 
                                                           em_min),
                                  csize=1000)

    # write log
    logfile = os.path.join(inferemdir, 'summary.log')
    if os.path.exists(inferemdir):
        f = open(logfile, 'a')
    else:
        f = open(logfile, 'w')
    
    f.write('\n*** start: Add Models to DB (Prior)***\n\n')
    f.write(arg_str)
    f.write('\n')
    f.write(summary_str)
    f.write('\n*** end: Add Models to DB (Prior)***\n')
    f.close()

    print summary_str
    
if __name__ == '__main__':
    # run
    main()

