# cmpy_bayes #

**Note:** Currently in heavy development.

A collection of Python scripts and a package of utility functions to make
Bayesian inference of unifilar hidden Markov models (uHMMs) and
epsilon-machines with Computation Mechanics Python, aka
[CMPy](http://cmpy.csc.ucdavis.edu/) easier.

## Dependencies ##

* [CMPy](http://cmpy.csc.ucdavis.edu/)

CMPy should be installed with the `--cython` switch to allow for use of the
enumeration library.  Also, note that CMPy has a quite a few dependencies that
are required to use these scripts.

## Installation ##

Get the code from github in the usual manner.  In a suitable directory clone
the repository using:

    git clone git@github.com:cstrelioff/cmpy_bayes.git

Install using:

    sudo python setup.py install

Or, for a local install, when administrator permissions are not available:

    python setup.py install --user

If you have issues with the scripts when using the `--user` switch, consult the
[local installation note](LOCALINSTALL.md).  The scripts are designed to be run
as system tools, allowing for execution from any directory without directly
calling python.

## Scripts ##

### cbayes_create_process_datafile.py ###

> Create a data file that works nicely with other cbayes scripts.  Can use
> any machine in `cmpy.machines` that has NO default parameters.  This means
> the process can be initiated using:
> 
> em = cmpy.machines.ProcessName()

### cbayes_enumerate_AddtoDB.py ###

