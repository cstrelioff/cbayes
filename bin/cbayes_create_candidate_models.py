#!/usr/bin/env python

"""cbayes_create_candidate_models.py


"""
from __future__ import division

import os
import argparse

import cmpy
import cmpy.inference.bayesianem as bayesem

from cbayes import create_dir
from cbayes import create_machine_file

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
    arg_list.append("-db     : Database directory "
            ">> {:s}\n".format(args.database_directory))
    arg_list.append("-a      : Alphabet size >> {:d}\n".format(args.alphabet_size))
    arg_list.append("-n      : Number of states "
            ">> {:s}\n".format(args.number_of_states))
    arg_list.append("--topological_eMs : "
           "topological eMs only? >> {:s}\n".format(str(args.topological_eMs)))
    arg_list.append("-nmax   : Max number of machines to load into RAM "
            ">> {:d}\n".format(args.nmax))
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
        """Create a `machines` file containing all candidate machine topologies
        to be used for inference."""
        )
    parser = argparse.ArgumentParser(description=desc_str)
    parser.add_argument('-db', '--database_directory',
            help = 'name of the database directory',
            type = create_dir,
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
    parser.add_argument('-nmax',
            help = 'maximum number of machines to load into RAM at once',
            type = int,
            default = 10000)
    parser.add_argument('-nprocs',
            help = 'number of simultaneous processes to run',
            type = int,
            default = 4)
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
    
    # create output file and write header
    machine_filename = os.path.join(args.database_directory, 'machines')
   
    # start the serious processing...
    summary_str = create_machine_file(machine_filename,
                                      args.alphabet_size,
                                      ns_list,
                                      em_min,
                                      args.nmax,
                                      args.nprocs)

    # write log
    logfile = os.path.join(args.database_directory, 'summary.log')
    if os.path.exists(logfile):
        f = open(logfile, 'a')
    else:
        f = open(logfile, 'w')
    
    f.write('\n*** start: Create set of candidate models for inference ***\n\n')
    f.write(arg_str)
    f.write('\n')
    f.write(summary_str)
    f.write('\n*** end: Create set of candidate models for inference ***\n')
    f.close()

    print summary_str
    
if __name__ == '__main__':
    main()

