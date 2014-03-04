#!/usr/bin/env python

"""cbayes_bash_enumerate_convergence.py

Create a bash script to run a sequence of `cbayes_enumerate_` scripts on a
specified data file.  In this case, the focus is on convergence of inference
by considering segments, using the single data file, of increasing length.

"""
import os
import shutil
import argparse
import numpy

from cbayes import check_dir_doesnot_exist
from cbayes import check_positive_int
from cbayes import check_positive_float
from cbayes import check_probability
from cbayes import create_dir
from cbayes import read_datafile

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
    parser.add_argument('--beta',
            help = 'list of beta values used for model comparison',
            type = str,
            required = True
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
    parser.add_argument('-nprocs',
            help = 'number of simultaneous processes to run',
            type = int,
            default = 4)

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
    arg_list.append("--beta  : Penalty size >> {}\n".format(args.beta))
    arg_list.append("-penalty : Type of penalty "
            ">> {:s}\n".format(args.penalty))
    arg_list.append("-sst : Subsample type >> "
                "{:s}\n".format(args.subsample_type))
    arg_list.append("-ns : Number of machines to sample "
            ">> {:d}\n".format(args.number_samples))
    arg_list.append("-nprocs : Number of simultaneous processes to run "
            ">> {:d}\n".format(args.nprocs))
    
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

    # process betas to consider
    beta_str = args.beta.split(',')
    beta_list = [float(b) for b in beta_str]

    ## write header
    header_list = []
    header_list.append("#!/bin/bash -l\n")
    header_list.append("## commands to execute\n\n")
    header_list.append("# keep track of total compute time\n")
    header_list.append("begin=\"$(date +%s)\"\n\n")
    f.write(''.join(header_list))

    ## calculate evidence terms for prior
    prior_list = []
    prior_list.append("##\n")
    prior_list.append("## PRIOR\n")
    prior_list.append("echo \">> Add PRIOR models to DB: `date`\"\n")
    # --single line
    prior_list.append("cbayes_enumerate_prior.py")
    prior_list.append(" -db {}".format(args.database_directory))
    prior_list.append(" -nprocs {}".format(args.nprocs))
    prior_list.append("\n")
    prior_list.append("echo\n")
    
    for beta in beta_list:
        ## calculate model prior probabilities
        prior_list.append("echo \">> Calculate PRIOR model "
                          "probabilities, beta= {} : `date`\"\n".format(beta))
        # --single line
        prior_list.append("cbayes_enumerate_probabilities.py")
        prior_list.append(" -db {}".format(args.database_directory))
        prior_list.append(" -idir inferEM_0-0") 
        prior_list.append(" --beta {}".format(beta)) 
        prior_list.append(" -p {}\n".format(args.penalty))
        prior_list.append("echo\n")
        ## sample machines from prior
        prior_list.append("echo \">> Sample PRIOR machines"
                ", beta = {} `date`\"\n".format(beta))
        # ---single line
        prior_list.append("cbayes_enumerate_sample.py")
        prior_list.append(" -db {}".format(args.database_directory))
        prior_list.append(" -idir inferEM_0-0")
        prior_list.append(" -mp probabilities_beta-{:.6f}".format(beta))
        prior_list.append("_penalty-{}".format(args.penalty))
        prior_list.append(" -ns {}".format(args.number_samples))
        prior_list.append(" --this_is_prior")
        prior_list.append(" -nprocs {}\n".format(args.nprocs))
        prior_list.append("echo\n")

    ## write to file
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
        posterior_list = []
        ## Add models to DB for this subsample
        posterior_list.append("##\n")
        posterior_list.append("echo \"SEGMENT : 0 -- {}\"\n".format(ssl))

        # --single line
        posterior_list.append("cbayes_enumerate_posterior.py")
        posterior_list.append(" -f {}".format(args.file))
        posterior_list.append(" -db {}".format(args.database_directory))
        posterior_list.append(" -sr 0,{}".format(ssl))
        posterior_list.append(" -nprocs {}".format(args.nprocs))
        posterior_list.append("\n")
        posterior_list.append("echo\n")

        for beta in beta_list:
            ## calculate model prior probabilities
            posterior_list.append("echo \">> Calculate POSTERIOR model "
                "probabilities, beta={} : `date`\"\n".format(beta))
            # --single line
            posterior_list.append("cbayes_enumerate_probabilities.py")
            posterior_list.append(" -db {}".format(args.database_directory))
            posterior_list.append(" -idir inferEM_0-{}".format(ssl)) 
            posterior_list.append(" --beta {}".format(beta)) 
            posterior_list.append(" -p {}\n".format(args.penalty))
            posterior_list.append("echo\n")
            ## sample machines from prior
            posterior_list.append("echo \">> Sample POSTERIOR machines,"
                    " beta={} : `date`\"\n".format(beta))
            # ---single line
            posterior_list.append("cbayes_enumerate_sample.py")
            posterior_list.append(" -db {}".format(args.database_directory))
            posterior_list.append(" -idir inferEM_0-{}".format(ssl))
            posterior_list.append(" -mp probabilities_beta-{:.6f}".format(beta))
            posterior_list.append("_penalty-{}".format(args.penalty))
            posterior_list.append(" -ns {}".format(args.number_samples))
            posterior_list.append(" -nprocs {}\n".format(args.nprocs))
            posterior_list.append("echo\n")

        ## write to file
        f.write(''.join(posterior_list))

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
    
    # check for existence of data and db/machines files
    check_dir_doesnot_exist(args.file)
    check_dir_doesnot_exist(os.path.join(args.database_directory, 'machines'))
    # copy data file to database_directory/datafile
    shutil.copy(args.file, os.path.join(args.database_directory, 'datafile'))

    # write scripts
    write_scripts(args)

if __name__ == '__main__':
    main()
