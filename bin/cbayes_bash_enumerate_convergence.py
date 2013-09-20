#!/usr/bin/env python

"""cbayes_bash_enumerate_convergence.py

Create a bash script to run a sequence of `cbayes_enumerate_` scripts on a
specified data file.  In this case, the focus is on convergence of inference
by considering subsegments, using the single data file, of increasing length.

"""
import os
import argparse
import numpy

from cmpy_bayes import check_positive_int
from cmpy_bayes import check_positive_float
from cmpy_bayes import check_probability
from cmpy_bayes import create_dir
from cmpy_bayes import read_datafile

def create_parser():
    """Create argparse instance for script.
    
    Returns
    -------
    args : agparse arg instance

    """
    desc_str = (
        """Run set of enumerative inferences on a single data file to look at
        convergence properties.
        
        """
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
            type = check_positive_int,
            default = 2,
            required = True)
    parser.add_argument('-n', '--number_of_states',
            help = 'a comma sep list of state numbers, eg 1,2,3',
            type = str,
            required = True)
    parser.add_argument('--beta',
            help = 'beta used for model comparison, penalty for num of states',
            type = check_positive_float,
            default = 4.0,
            required = False
            )
    parser.add_argument('-p', '--penalty',
            help = 'type of penalty to apply',
            type = str,
            choices = ['num_states', 'num_edges'],
            default = 'num_states',
            required = True)
    parser.add_argument('-sst', '--subsample_type',
            help = ("subsample type, for monitoring convergence"),
            type = str,
            default = 'powers_of_2',
            choices = ['powers_of_2', 'powers_of_10'],
            required = False
            )
    parser.add_argument('-ns', '--number_samples',
            help = 'number of machines to sample',
            type = check_positive_int,
            required = True
            )
    parser.add_argument('--topological_eMs',
            help = 'include only topological eM structures',
            action = 'store_true',
            required = False
            )

    # do the parsing
    args = parser.parse_args()

    return args

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
    arg_list.append("--beta  : Penalty size >> {:f}\n".format(args.beta))
    arg_list.append("-penalty : Type of penalty "
            ">> {:s}\n".format(args.penalty))
    arg_list.append("--topological_eMs : "
           "topological eMs only? >> {:s}\n".format(str(args.topological_eMs)))
    arg_list.append("-sst : Subsample type >> "
                "{:s}\n".format(args.subsample_type))
    arg_list.append("-ns : Number of machines to sample "
            ">> {:d}\n".format(args.number_samples))
    
    arg_str = ''.join(arg_list)

    return arg_str

def write_scripts(args):
    """Write the slurm script.
    
    Parameters
    ----------
    args : results of argparse
        Command line arguments.
    outdir : str
        Output (database) directory.

    """
    # create filename and open
    jobname ="enumerate_bash_converge_{:s}".format(args.file.split('.')[0])
    fname = ''.join([jobname, ".sh"])
    f = open(fname, 'w')

    ## write header
    header_list = []
    header_list.append("#!/bin/bash -l\n")
    header_list.append("## commands to execute\n\n")
    header_list.append("# keep track of total compute time\n")
    header_list.append("begin=\"$(date +%s)\"\n\n")
    f.write(''.join(header_list))

    ## add prior to script - analyze machines, prior
    prior_list = []
    prior_list.append("##\n")
    prior_list.append("## PRIOR\n")
    prior_list.append("echo \">> Add PRIOR models to DB: `date`\"\n")
    # single line
    prior_list.append("cbayes_enumerate_PriorAddToDB.py")
    prior_list.append(" -db {}".format(args.database_directory))
    prior_list.append(" -a {} ".format(args.alphabet_size))
    prior_list.append(" -n {}".format(args.number_of_states))
    if args.topological_eMs:
        prior_list.append(" --topological_eMs\n")
    else:
        prior_list.append("\n")
    prior_list.append("echo\n")
    prior_list.append("echo \">> Calculate PRIOR model "
            "probabilities: `date`\"\n")
    
    # single line - calculate probability for machines using prior
    prior_list.append("cbayes_enumerate_CalcProbs.py")
    prior_list.append(" -db {}".format(args.database_directory))
    prior_list.append(" -idir inferEM_0-0") 
    prior_list.append(" --beta {}".format(args.beta)) 
    prior_list.append(" -p {}\n".format(args.penalty))
    prior_list.append("echo\n")
    prior_list.append("echo \">> Sample PRIOR machines: `date`\"\n")
    
    # single line - sample machines using prior
    prior_list.append("cbayes_enumerate_Sample.py")
    prior_list.append(" -db {}".format(args.database_directory))
    prior_list.append(" -idir inferEM_0-0")
    prior_list.append(" -mp modelprobs_beta-{:.6f}".format(args.beta))
    prior_list.append("_penalty-{}.pickle".format(args.penalty))
    prior_list.append(" -ns {}\n".format(args.number_samples))
    prior_list.append("echo\n")
    prior_list.append("echo \">> Process sampled PRIOR machines: `date`\"\n")

    # single line - process sampled machines
    prior_list.append("cbayes_enumerate_ProcessSamples.py")
    prior_list.append(" -db {}".format(args.database_directory))
    prior_list.append(" -sdir samples_0-0")
    prior_list.append("_beta-{:.6f}".format(args.beta))
    prior_list.append("_penalty-{}\n".format(args.penalty))
    prior_list.append("echo\n")
    f.write(''.join(prior_list))
    
    # read data to get data length
    data = read_datafile(args.file)
    data_len = len(data)
    del data
    
    ## Add subsections of data for convergence analysis
    num_steps = 0
    if args.subsample_type == 'powers_of_2':
        num_steps = numpy.floor(numpy.log2(data_len))
        data_length_list = [2**n for n in range(0, int(num_steps)+1)]
    elif args.subsample_type == 'powers_of_10':
        num_steps = numpy.floor(numpy.log10(data_len))
        data_length_list = [10**n for n in range(0, int(num_steps)+1)]

    ## iterate through subsample length, add to script
    for ssl in data_length_list:
        ssl_list = []
        ## Add models to DB for this subsample length
        ssl_list.append("##\n")
        ssl_list.append("## SUBSAMPLE, LENGTH: {}\n".format(ssl))
        ssl_list.append("echo \">> Add models, subsample")
        ssl_list.append(" length {}, to DB: `date`\"\n".format(ssl))
        # single line
        ssl_list.append("cbayes_enumerate_AddToDB.py")
        ssl_list.append(" -f {}".format(args.file))
        ssl_list.append(" -db {}".format(args.database_directory))
        ssl_list.append(" -a {} ".format(args.alphabet_size))
        ssl_list.append(" -n {}".format(args.number_of_states))
        ssl_list.append(" -sr 0,{}".format(ssl))
        if args.topological_eMs:
            ssl_list.append(" --topological_eMs\n")
        else:
            ssl_list.append("\n")

        ## calculate model probabilities for the subsample length
        ssl_list.append("echo\n")
        ssl_list.append("echo \">> Calculate model "
                "probabilities: `date`\"\n")
        # single line
        ssl_list.append("cbayes_enumerate_CalcProbs.py")
        ssl_list.append(" -db {}".format(args.database_directory))
        ssl_list.append(" -idir inferEM_0-{}".format(ssl)) 
        ssl_list.append(" --beta {}".format(args.beta)) 
        ssl_list.append(" -p {}\n".format(args.penalty))

        ## sample machines for this subsample length
        ssl_list.append("echo\n")
        ssl_list.append("echo \">> Sample machines: `date`\"\n")
        
        # single line -- sample machines
        ssl_list.append("cbayes_enumerate_Sample.py")
        ssl_list.append(" -db {}".format(args.database_directory))
        ssl_list.append(" -idir inferEM_0-{}".format(ssl))
        ssl_list.append(" -mp modelprobs_beta-{:.6f}".format(args.beta))
        ssl_list.append("_penalty-{}.pickle".format(args.penalty))
        ssl_list.append(" -ns {}\n".format(args.number_samples))
        
        ## process sampled machines
        ssl_list.append("echo\n")
        ssl_list.append("echo \">> Process sampled machines: `date`\"\n")
        
        # single line - process sampled machines
        ssl_list.append("cbayes_enumerate_ProcessSamples.py")
        ssl_list.append(" -db {}".format(args.database_directory))
        ssl_list.append(" -sdir samples_0-{}".format(ssl))
        ssl_list.append("_beta-{:.6f}".format(args.beta))
        ssl_list.append("_penalty-{}\n".format(args.penalty))
        ssl_list.append("echo\n")

        # write to file for this subsample length
        f.write(''.join(ssl_list))

    # calculate total compute time
    f.write("# calculate total compute time\n")
    f.write("end=\"$(date +%s)\"\n")
    f.write("diff=\"$(($end-$begin))\"\n")
    f.write("printf \"Total Compute Time: %02d:%02d:%02d:%02d\"" 
            " \"$((diff/86400))\" \"$((diff/3600%24))\"" 
            " \"$((diff/60%60))\" \"$((diff%60))\"\n")
    f.write("echo\n")
    f.write("echo\n")
    f.close()

def main():
    """Create bash script for given data file and parameters."""

    # get command line args
    args=create_parser()

    # reports args
    summary_str = report_args(args)
    print summary_str

    # write scripts
    write_scripts(args)

if __name__ == '__main__':
    main()
