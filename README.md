# cmpy_bayes modules #

**Note:** Currently in heavy development.

A collection of Python scripts and a package of utility functions to make
Bayesian inference of unifilar hidden Markov models (uHMMs) and
epsilon-machines with Computation Mechanics Python, aka
[CMPy](http://cmpy.csc.ucdavis.edu/) easier.

## Scripts ##

### cbayes_create_process_datafile.py ###

    Create a data file that works nicely with other cbayes scripts.  Can use
    any machine in `cmpy.machines` that has NO default parameters.  This means
    the process can be initiated using:

    em = cmpy.machines.ProcessName()

### cbayes_enumerate_AddtoDB.py ###

## Dependencies ##

CMPy -- Developed at the 

## Install ##

Get the code from github in the usual manner.  In a suitable directory clone
the respository using:

    git clone git@github.com:cstrelioff/research.git

Install using:

    sudo python setup.py install

Or, for local install, when admin permission are not available:

    python setup.py install --user

This command will install to .local/bin in the user home directory on Ubuntu.
This location is not in the system path by default, so we add this location to
the path.  In the .bashrc file (create a file in your home directory if it does
not already exist) add the following lines:

    # include .local/bin for local python scripts
    export PATH=~/.local/bin/:$PATH

To make this definition active for the current session you have to `source` the
.bashrc file using:

    source .bashrc

Turns out that an interactive login does not source `.bashrc`.  So, create a
`.profile` file in your server home directory with the following contents:

    # ~/.profile: executed by the command interpreter for login shells.
    # This file is not read by bash, if 
    #    ~/.bash_profile 
    # or 
    #    ~/.bash_login
    # exist
    
    # if running bash
    if [ -n "$BASH_VERSION" ]; then
        # include .bashrc if it exists
        if [ -f "$HOME/.bashrc" ]; then
            . "$HOME/.bashrc"
        fi
    fi
    
    # set PATH so it includes user's private bin if it exists
    if [ -d "$HOME/bin" ] ; then
        PATH="$HOME/bin:$PATH"
    fi

## Scripts ##

Modules/scripts included in this repository are:

### cbayes\_gabayes.py ###

This evolves the epsilon-machine structure that best matches a provided data
file.  Bayesian model comparison is applied, taking the place of fitness in a
normal (GA) genetic algorithm.

Get help using

    cbayes_gabayes.py --help





