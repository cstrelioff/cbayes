#!/usr/bin/env python


import os
import argparse
import numpy

from cbayes import check_positive_int
from cbayes import check_positive_float
from cbayes import check_probability
from cbayes import create_dir

def create_parser():
    """Create argparse instance for script.
    
    Returns
    -------
    args : agparse arg instance

    """
    desc_str = (
        """Run lots of ga_bayes jobs on slurm."""
        )

    parser = argparse.ArgumentParser(description=desc_str)
    parser.add_argument('-f','--file', 
            help = 'input (data) file name', 
            type = str,
            required = True)
    parser.add_argument('-r', '--runs_number',
            help = 'number of runs for given data file',
            type = check_positive_int,
            required = True)
    parser.add_argument('-a', '--alphabet_size',
            help = 'number of letters in alphabet',
            type = int,
            default = 2,
            required = True)
    parser.add_argument('-g', '--generations',
            help = 'generations to run GA',
            type = check_positive_int,
            required = True)
    parser.add_argument('-d','--directory', 
            help = 'output directory for all ga_bayes runs', 
            type = str,
            required = True)
    parser.add_argument('-beta',
            help = 'beta used for model comparison, penalty for num of states',
            type = check_positive_float,
            default = 4.0,
            required = False
            )
    parser.add_argument('-n', '--population_size',
            help = 'number of individuals in evolving population',
            type = int,
            default = 500,
            required = False
            )
    parser.add_argument('--population_selection',
            help = 'selection using population effects?',
            action = 'store_true',
            default = False
            )
    parser.add_argument('-pi', '--initial_population',
            help = 'makeup of the inital population',
            type = str,
            default = 'single-state',
            choices = ['single-state', 'library'],
            required = False
            )
    parser.add_argument('-mu_p', '--mutation_point',
            help = 'per-site, per-generation mutation rate',
            type = check_probability,
            default = 0.1,
            required = False)
    parser.add_argument('-mu_i','--mutation_insert',
            help = 'per-generation >insert< mutation rate',
            type = check_probability,
            default = 0.1,
            required = False)
    parser.add_argument('-mu_d', '--mutation_delete',
            help = 'per-generation >deletion< mutation rate',
            type = check_probability,
            default = 0.1,
            required = False)
    parser.add_argument('--variable_insert',
            help = 'allow variable length insert mutations?',
            action = 'store_true',
            default = False)
    parser.add_argument('--generation_summary',
            help = 'generation summary?',
            action = 'store_true',
            default = False)

    # do the parsing
    args = parser.parse_args()

    return args

def report_args(args):
    """Report the requested settings.
    
    Returns
    -------
    summary_str : str
        A summary of run settings.

    """

    summary_list = []
    summary_list.append("SETTINGS:\n")
    summary_list.append("-f    : Data file >> {:s}\n".format(args.file))
    summary_list.append("-r    : Number of runs "
            ">> {:d}\n\n".format(args.runs_number))
    summary_list.append("-a    : Alphabet size "
            ">> {:d}\n".format(args.alphabet_size))
    summary_list.append("-g    : Generations to evolve "
            ">> {:d}\n".format(args.generations))
    summary_list.append("-d    : Output directory file "
            ">>  {:s}\n".format(args.directory))
    summary_list.append("-beta : for model size penalty "
            ">>  {:f}\n".format(args.beta))
    summary_list.append("-n    : Population size "
            ">> {:d}\n".format(args.population_size)) 
    summary_list.append("--population_selection: Population selection? "
            ">> {:s}\n".format(str(args.population_selection)))
    summary_list.append("-pi   : Initial population "
            ">> {:s}\n".format(args.initial_population))
    summary_list.append("-mu_p : mutation, pt "
            ">> {:f}\n".format(args.mutation_point))
    summary_list.append("-mu_i : mutation, insert "
            ">> {:f}\n".format(args.mutation_insert))
    summary_list.append("-mu_d : mutation, delete "
            ">> {:f}\n".format(args.mutation_delete))
    summary_list.append("--variable_insert: variable length? "
            ">> {:s}\n".format(str(args.variable_insert)))
    summary_list.append("--generation_summary: summary file for"
        " each generation? >> {:s}\n".format(str(args.generation_summary)))
    summary_list.append("\n")
    
    summary_str = ''.join(summary_list)

    return summary_str

def write_scripts(args, outdir):
    """Write the slurm script."""

    # change to outdir, remembering cwd
    cwd=os.getcwd()
    os.chdir(outdir)

    # create rng for random seeds
    rng = numpy.random

    for n in range(args.runs_number):
        # setup script string for format below
        script = ("#!/bin/bash -l\n"
        "# NOTE the -l flag!\n\n"
        "# Name of the job\n"
        "#SBATCH -J {jobname:s}\n\n"
        "# Standard out and Standard Error output files\n"
        "#SBATCH -o {jobname:s}-%j.output\n"
        "#SBATCH -e {jobname:s}-%j.output\n\n"
        "## commands to execute\n\n"
        "echo\n"
        "echo start: `date`\n"
        "srun -l research_gabayes.py "
        "-f ../{file:s} "
        "-a {alphabet_size:d} "
        "-g {generations:d} "
        "-s {seed:d} "
        "-d {directory:s} "
        "-beta {beta:f} "
        "-n {population_size:d} "
        "-pi {initial_population:s} "
        "-mu_p {mutation_point:f} "
        "-mu_i {mutation_insert:f} "
        "-mu_d {mutation_delete:f} "
        "{flags:s}"
        "\n"
        "echo end: `date`\n"
        "echo\n"
        )
        
        # create file
        fstart = "r{:03d}-".format(n) + args.file.split('.')[0]
        f = open(fstart + '.sh', 'w')

        # parse flags
        flags=[]
        if args.variable_insert:
            flags.append("--variable_insert")
        if args.population_selection:
            flags.append("--population_selection")
        if args.generation_summary:
            flags.append("--generation_summary")

        # create string
        flag_str = ' '.join(flags)

        # populate script
        script = script.format(jobname=fstart,
                file=args.file,
                alphabet_size=args.alphabet_size,
                generations=args.generations,
                directory=fstart,
                seed=rng.randint(0,1000000),
                beta=args.beta,
                population_size=args.population_size,
                initial_population=args.initial_population,
                mutation_point=args.mutation_point,
                mutation_insert=args.mutation_insert,
                mutation_delete=args.mutation_delete,
                flags=flag_str
                )

        f.write(script)
        f.close()

    # return to cwd
    os.chdir(cwd)

def main():
    """Create slurm scripts for given data file and parameters."""

    # get command line args
    args=create_parser()

    # reports args
    summary_str = report_args(args)
    print summary_str

    # create ouput directory
    outdir = create_dir(args.directory)

    # write scripts
    write_scripts(args, outdir)

if __name__ == '__main__':
    main()
