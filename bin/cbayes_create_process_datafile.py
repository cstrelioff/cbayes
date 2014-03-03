#!/usr/bin/env python

"""cbayes_create_process_datafile.py

Create a data file that works nicely with other cbayes scripts.  Can use
any machine in `cmpy.machines` that **has default parameters**.  This means
the process can be instantiated using:

em = cmpy.machines.ProcessName()

"""

from __future__ import division

import os
import argparse

from cbayes import write_data
from cbayes import write_em_pickle

try:
    import cmpy
    import cmpy.inference.bayesianem as bayesem
    import cmpy.orderlygen.pyicdfa as pyicdfa
except:
    raise Exception("CMPy must be installed with --cython switch!")

# exception
class ProcessException(Exception):
    pass

def create_process(args):
    """Check that specified process makes sense and create eM.
    
    Parameter
    ---------
    args : object
        Arguments from command line argparse.
    
    Returns
    -------
    eM : RecurrentEpsilonMachine
        A recurrent eM for the desired process.

    Exception
    ---------
        Exception raised if process requested does not make sense.

    """
    process = args.process

    # create list of valid machines
    valid_machines = []
    valid_types = [cmpy.machines.MealyHMM, 
            cmpy.machines.RecurrentEpsilonMachine]

    for em in dir(cmpy.machines):
        if em[0].isupper():
            try:
                m_str = 'cmpy.machines.' + em +'()' 
                eval(m_str)
                mtype = type(eval(m_str))
                if mtype in valid_types:
                    valid_machines.append(em)
            except:
                pass

    # remove MealyHMM, RecurrentEpsilonMachine
    valid_machines.remove('MealyHMM')
    valid_machines.remove('RecurrentEpsilonMachine')

    # if in valid_machine, try to create instance
    if process in valid_machines:
        eM = eval('cmpy.machines.' + process + '()')
    else: 
        error_msg = ("\n\nProcess {} not valid. Try:\n\n{}\n".format(process,
                                                              valid_machines))
        raise ProcessException(error_msg)

    return eM

def create_parser():
    """Create argparse instance for script.
    
    Returns
    -------
    args : agparse arg instance

    """
    desc_str = ("""\nCreate a datafile for the specified process.\n"""
        )

    parser = argparse.ArgumentParser(description=desc_str)
    parser.add_argument('-f','--file', 
            help = 'output (data) file name', 
            type = str,
            required = True)
    parser.add_argument('-l', '--length',
            help = 'length of the time series',
            type = int,
            required = True)
    parser.add_argument('-p', '--process',
            help = ("process from cmpy.machines, eg Even"),
            type = str,
            required = True
            )
    parser.add_argument('-d','--directory', 
            help = 'output directory name', 
            type = str,
            default = '.',
            required = False)

    # do the parsing
    args = parser.parse_args()

    return args

def report_args(args):
    """Report the requested settings."""

    print ("SETTINGS:\n")
    print ("-f    : Output data file >> {:s}".format(args.file))
    print ("-l    : Length of data series >> {:d}".format(args.length))
    print ("-p    : Process >> {:s}".format(args.process))
    print ("-d    : Ouput diretory >> {:s}".format(args.directory))
    print ("\n")

def main():
    """Create a datafile for specified process."""
    # get command line args
    args = create_parser()

    # report args
    report_args(args)

    # check and create instance of process, if possible
    eM = create_process(args)

    # write data
    write_data(args.directory, args.file, args.process, args.length, eM)

    # write machine to pickle
    write_em_pickle(args.file, eM)

if __name__ == '__main__':
    # run
    main()
