#!/usr/bin/env python

import os
import argparse

import cmpy
import cmpy.inference.evolveem as evolveem

from research import check_dir_exists
from research import check_positive_float
from research import check_positive_int
from research import check_probability
from research import read_datafile

def evolve(ev_pop, generations, output_dir, viable_only=True, verbose=False,
        generation_summary=False):
    """Run an instance of the EvolveEM class with additional paramters.
    
    Paramerters
    -----------
    ev_pop : EvolveEM
        An instance of the EvolveEM class.
    generations : int
        The number of generations to run the evolutionary algorithm.
    output_dir : str
        Name of the output directory for all files and plots created.
    viable_only : bool [default True]
        Only write viable genomes to generation summaries, if the files are
        written.
    verbose : bool [default False]
        Provide lots per-generation summary to stdout.
    generation_summary : bool [default False]
        If true, write a summary file for each generation.

    """

    ## generation 0
    if generation_summary:
        ev_pop.write_generation_summary(dirname=output_dir,
                                        viable_only=viable_only)
    ev_pop.draw_most_probable(dirname=output_dir)
    ev_pop.pickle_most_probable(dirname=output_dir)
    if verbose:
        print ev_pop.summary_string()

    ## generations 1-
    for n in range(1,generations+1):
        ev_pop.next_generation()
        if generation_summary:
            ev_pop.write_generation_summary(dirname=output_dir, 
                                            viable_only=viable_only)
        ev_pop.draw_most_probable(dirname=output_dir)
        ev_pop.pickle_most_probable(dirname=output_dir)
        if verbose:
            print ev_pop.summary_string()

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
    summary_list.append("-a    : Alphabet size "
            ">> {:d}\n".format(args.alphabet_size))
    summary_list.append("-s    : Seed for RNG >> {:d}\n".format(args.seed))
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

def create_parser():
    """Create argparse instance for script.
    
    Returns
    -------
    args : agparse arg instance

    """
    desc_str = (
        """Evolve epsilon-machine structure for supplied data."""
        )

    parser = argparse.ArgumentParser(description=desc_str)
    parser.add_argument('-f','--file', 
            help = 'input (data) file name', 
            type = str,
            required = True)
    parser.add_argument('-a', '--alphabet_size',
            help = 'number of letters in alphabet',
            type = check_positive_int,
            default = 2,
            required = True)
    parser.add_argument('-g', '--generations',
            help = 'generations to run GA',
            type = check_positive_int,
            required = True)
    parser.add_argument('-s', '--seed',
            help = 'seed for random number generator',
            type = check_positive_int,
            required = True)
    parser.add_argument('-d','--directory', 
            help = 'output directory name', 
            type = str,
            default = 'EvolveEMOut',
            required = True)
    parser.add_argument('-beta',
            help = 'beta used for model comparison, penalty for num of states',
            type = check_positive_float,
            default = 4.0,
            required = False
            )
    parser.add_argument('-n', '--population_size',
            help = 'number of individuals in evolving population',
            type = check_positive_int,
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

def initialize_rng(seed):
    """Initial rng with passed seed.
    
    Parameters
    ----------
    seed : int
        Seed for random number generator.
        
    Returns
    -------
    rng : RNG
        A random number generator.
        
    """
    import cmpy.math
    rng = cmpy.math.prng
    rng.seed(seed)

    return rng

def initialize_ga_pop(data, prng, args):
    """Initializ instance of EvolvePopulationEM.

    Parameters
    ----------
    data : list
        A data series as a list.
    prng : RNG
        A rng, proabably created using initialize_rng function.
    args : argparse args
        Reulst of command line argparse.

    Returns
    -------
        ev_pop : instance of EvolvePopulationEM

    """
    ev_pop = evolveem.EvolvePopulationEM(data, 
                                    args.alphabet_size,
                                    args.population_size, 
                                    initial_population=args.initial_population,
                                    mutation_pt=args.mutation_point,
                                    mutation_insert=args.mutation_insert, 
                                    mutation_delete=args.mutation_delete,
                                    var_length_insert=args.variable_insert,
                                    pop_selection=args.population_selection, 
                                    beta=args.beta, 
                                    summary_file=True,
                                    dir_name=args.directory,
                                    prng=prng)
    return ev_pop


def main():
    """Run the evolutionary inference algorithm with passed data,
    parameters.
    
    """
    # parse command line
    args = create_parser()

    # get command line args and report settings
    summary_str = report_args(args)
    print summary_str

    # check directory
    check_dir_exists(args.directory)

    # get rng
    prng = initialize_rng(args.seed)

    # read file
    data = read_datafile(args.file)

    # initalize
    ev_pop = initialize_ga_pop(data, prng, args)

    # run
    evolve(ev_pop, args.generations, args.directory, viable_only=True,
            verbose=False, generation_summary=args.generation_summary)

if __name__ == '__main__':
    main()

