#!/usr/bin/env python

"""cbayes_enumerate_probabilities.py

This script processes the output of cbayes_enumerate_prior.py or
cbayes_enumerate_posterior.py scripts to calculate prior, or posterior,
probabilities for each machine/model topolgy.  This script must be run before
samples can be generated from the prior or posterior.

"""
from __future__ import division

import os
import argparse

import cmpy
import cmpy.inference.bayesianem as bayesem

from cbayes import check_positive_float
from cbayes import read_datafile
from cbayes import calc_probs_beta_db

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
    arg_list.append("-idir : InferEM  sub-directory "
            ">> {:s}\n".format(args.inferem_directory))
    arg_list.append("-beta  : Penalty size >> {:f}\n".format(args.beta))
    arg_list.append("-penalty : Type of penalty "
            ">> {:s}\n".format(args.penalty))
    
    arg_str = ''.join(arg_list)

    return arg_str

def create_parser():
    """Create argparse instance for script.
    
    Returns
    -------
    args : agparse arg instance

    """
    desc_str = (
        """Calculate the model probabilities for a directory of pickled
        InferEM instances using the specified penalty type and strength."""
        )
    parser = argparse.ArgumentParser(description=desc_str)
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
    parser.add_argument('--beta',
            help = 'strength of exponential penalty',
            type = check_positive_float,
            required = True)
    parser.add_argument('-p', '--penalty',
            help = 'type of penalty to apply',
            type = str,
            choices = ['num_states', 'num_edges'],
            default = 'num_states',
            required = True)
    
    # do the parsing
    args = parser.parse_args()

    return args

def main():
    """Calculate the prior, or posterior, probabilities for the specified
    inferEM_*-* directory.
    
    """
    # parse command line
    args = create_parser()

    # get command line args and report settings
    arg_str = report_args(args)
    print arg_str

    # generate model probabilities for models in DB 
    summary_str, model_probs = calc_probs_beta_db(args.database_directory,
                                    args.inferem_directory,
                                    args.beta,
                                    args.penalty)

    # write log
    logfile = os.path.join(args.database_directory, 'summary.log')
    if os.path.exists(logfile):
        f = open(logfile, 'a')
    else:
        f = open(logfile, 'w')
     
    f.write('\n*** start: calculate model probabilities***\n\n')
    f.write(arg_str)
    f.write('\n')
    f.write(summary_str)
    f.write('\n*** end: calculate model probabilities***\n')
    f.close()
    
    print summary_str
    print "Top 10 machines"
    print "\n{:>15s} {:>15s} {:>15s}".format('Model', 'Prob',
                                          'Cumm Prob')
    cumm_pr = 0.
    for eM in sorted(model_probs,key=model_probs.get, reverse=True)[0:10]:
        cumm_pr += model_probs[eM]
        print "{:>15s} {:>15e} {:>15e}".format(eM, model_probs[eM], cumm_pr)
    
if __name__ == '__main__':
    main()

