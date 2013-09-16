#!/usr/bin/env python

"""cbayes_create_process_datafile.py

Christoper C. Strelioff
chris.strelioff@gmail.com

Create a data file that works nicely with other cbayes scripts.  Can use
any machine in `cmpy.machines` that **has default parameters**.  This means
the process can be instantiated using:

em = cmpy.machines.ProcessName()

"""

from __future__ import division

import os
import argparse

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

def write_data(args, eM):
    """Write the data file."""
    
    cwd = os.getcwd()
    out_dir = cwd
    if args.directory != '.':
        out_dir = os.path.join(cwd, args.directory)
        if os.path.exists(out_dir):
            os.chdir(out_dir)
        else:
            raise ProcessException("Output directory does not exist.")

    f = open(args.file, 'w')

    # header
    f.write("# SETTINGS:\n")
    f.write("# -f    : Output data file >> {:s}\n".format(args.file))
    f.write("# -l    : Length of data series >> {:d}\n".format(args.length))
    f.write("# -p    : Process >> {:s}\n".format(args.process))
    f.write("# -d    : Ouput diretory >> {:s}\n".format(args.directory))
    f.write('#\n')

    # eM details
    f.write('# eM \n')
    f.write("# name: {:s}\n".format(str(eM)))

    trans_str = eM.to_string()
    trans_list = trans_str.split(';')
    for element in trans_list:
        element = element.strip()
        f.write('# ' + element + '\n')

    # now, data
    eM.randomize_current_node()
    sN = eM.get_current_node()
    f.write("# start node: {:s}\n".format(str(sN)))
    f.write('#\n')
    
    # create the data
    data = eM.symbols(args.length)
    data_str = ','.join(data)

    # find the evidence for generated data
    evidence_available = True
    try:
        posterior = bayesem.InferEM(eM, data)
        log_evidence = posterior.log_evidence()
        for beta in [0.0, 2.0, 4.0]:
            curr_evi = log_evidence - beta*len(eM.nodes())
            tmp= ("# log P(D|Mi)-beta*|Mi| "
                "(beta={b:f}) = {ce:.4f}\n"
                )
            f.write(tmp.format(b=beta, ce=curr_evi))
        f.write("#\n")
    except:
        evidence_available = False
        #f.write('# log P(D|Mi) not found -- non-unifilar source??\n')

    # compare with IID
    iid = cmpy.machines.IID(len(eM.alphabet()))
    iid_posterior = bayesem.InferEM(iid, data)
    iid_evidence = iid_posterior.log_evidence()
    f.write("# log P(D|IID) (beta=0.0) = {:.4f}\n".format(iid_evidence))
    f.write("#\n")

    # write data and close file
    f.write(data_str)
    f.close()

    # include an inference log file
    if evidence_available:
        flog = open(args.file.split('.')[0] + '.log', 'w')
        flog.write("# Inference information for {:s}\n".format(args.file))
        flog.write("Evi_True_beta_0 Evi_True_beta_2 ")
        flog.write("Evi_True_beta_4 Evi_IID_beta_0\n")
        flog.write("{:.2f} ".format(log_evidence))
        flog.write("{:.2f} ".format(log_evidence - 2.*len(eM.nodes())))
        flog.write("{:.2f} ".format(log_evidence - 4.*len(eM.nodes())))
        flog.write("{:.2f}\n".format(iid_evidence))
        flog.close()
    else:
        flog = open(args.file.split('.')[0] + '.log', 'w')
        flog.write("# Inference information for {:s}\n".format(args.file))
        flog.write("Evi_IID_beta_0\n")
        flog.write("{:.2f}\n".format(iid_evidence))
        flog.close()

    # change back to entry directory
    os.chdir(cwd)

def write_em_pickle(args, eM):
    """Write a pickle file of the generated process for later use.
    
    Parmeters
    ---------
    args : result of argarse
        Results of command line argument parsing.
    eM : RecurrentEpsilonMachine, MealyHMM
        Machine to pickle.

    """
    try:
        import cPickle as pickle
    except:
        import pickle
    
    try:
        filename = args.file.split('.')[0]
    except:
        filename = args.file
    f = open(filename + '.pickle', 'w')
    pickle.dump(eM, f)
    f.close()


def main():
    """Create a datafile for specified process."""
    # get command line args
    args = create_parser()

    # report args
    report_args(args)

    # check and create instance of process, if possible
    eM = create_process(args)

    # write data
    write_data(args, eM)

    # write machine to pickle
    write_em_pickle(args, eM)

if __name__ == '__main__':
    # run
    main()
