# cbayes #

**Note:** Currently in heavy development.

A collection of Python scripts to make Bayesian inference of unifilar hidden
Markov models (uHMMs) and epsilon-machines with Computational Mechanics Python,
aka [CMPy](http://cmpy.csc.ucdavis.edu/) easier.  The goal is to have an
easy-to-use pipeline of scripts that maximize the use of all compute power on a
single machine-- this is done using
[multiprocessing](http://docs.python.org/2/library/multiprocessing.html) to
launch many processes whenever possible.

## Dependencies ##

* [CMPy](http://cmpy.csc.ucdavis.edu/)

CMPy should be installed with the `--cython` switch to allow for use of the
enumeration library.  Also, note that CMPy has a quite a few
[dependencies](http://cmpy.csc.ucdavis.edu/installation.html) that
are required to use these scripts.

## Installation ##

Get the code from github in the usual manner.  In a suitable directory clone
the repository using:

    git clone git@github.com:cstrelioff/cbayes.git

Install using:

    sudo python setup.py install

Or, for a local install, when administrator permissions are not available:

    python setup.py install --user

If you have issues with the scripts when using the `--user` switch, consult the
[local installation note](LOCALINSTALL.md).  The scripts are designed to be run
as system tools, allowing for execution from any directory without directly
calling python.

## Scripts ##

All of these scripts should install as system scripts, making use from the
command line easy.  To get help for a script, type

    scriptname --help

### High Level Scripts ###

The scripts produce a *master script* that calls each of the base scripts (see
below) repeatedly to accomplish some larger inference goal.  Broadly speaking
there are two types of analyses:

* **Convergence** -- a single data series is repeatedly analyzed, starting with
  a small segment, and increasing the segment length until the full data set is
  used.  This analysis is helpful for looking at convergence of inference with
  increasing data size.

* **Overlap** -- a single data series is repeatedly analyzed, using reasonable
  sized overlapping chunks, to test for stationary behavior.  The same model
  topologies should probable for each segment, reflecting in-class and
  stationary behavior.

The **bash** scripts are designed to run on any Linux/Mac machine and the
**slurm** scripts can be used on a cluster running the
[slurm](https://computing.llnl.gov/linux/slurm/) resource manager.

#### cbayes_bash_enumerate_convergence.py ###

Create a **bash** script to run a sequence of `cbayes_enumerate_` scripts on a
specified data file.  In this case, the focus is on convergence of inference
by considering segments, using the single data file, of increasing length.

#### cbayes_bash_enumerate_overlap_analysis.py ####

Creates a **bash** script to run a sequence of `cbayes_enumerate_` scripts on a
specified data file.  In this case, the focus is on analysis of single data
file as well as considering overlapping subsegments to test for stationary
behavior.  In this context, stationary means that analysis of segments of the
total data series returns the same model(s) at different points -- reflecting a
single model topology, or set of topologies, is appropriate for describing the
data series.

#### cbayes_slurm_enumerate_convergence.py ####

Create a **slurm** script to run a sequence of `cbayes_enumerate_` scripts on a
specified data file.  In this case, the focus is on convergence of inference
by considering subsegments, using the single data file, of increasing length.

#### cbayes_slurm_enumerate_overlap_analysis.py ####

Create a **slurm** script to run a sequence of `cbayes_enumerate_` scripts on a
specified data file.  In this case, the focus is on analysis of single data
file as well as considering overlapping segments to test for stationary
behavior.  In this context, stationary means that analysis of segments of the
total data series return the same model at different points -- reflecting a
single, static model topology is appropriate for the complete data series.

### Base Scripts ###

These base scripts accomplish very specific tasks.  While they can be used alone,
it is often useful to employ the high level scripts (see above) to generate a
bash or slurm master script for generation of a whole research protocol.  These
master scripts use base scripts described below.

#### cbayes_create_candidate_models.py ####

This script is used to create a file with all candidate machine topologies
that the user wants to consider.  In addition, a database directory -- really
just a normal directory -- will be created for later use by other scripts.


#### cbayes_create_process_datafile.py ####

Create a data file that works nicely with other cbayes scripts.  Can use
any machine in `cmpy.machines` that **has default parameters**.  This means
the process can be instantiated using:

    em = cmpy.machines.ProcessName()

**Note:** If you have your own data, the output of this script would be a good
reference to see how the data file is formatted.


#### cbayes_enumerate_prior.py ####

This script focuses on the prior over models, calculating the prior log evidence
for all model topologies in the `machines` file.

#### cbayes_enumerate_posterior.py ####

This script focuses on the posterior over models, calculating the posterior log
evidence for all model topologies in the `machines` file using the specified
data file and data range.


#### cbayes_enumerate_probabilities.py ####

This script processes the output of `cbayes_enumerate_prior.py` or
`cbayes_enumerate_posterior.py` scripts to calculate prior, or posterior,
probabilities for each machine/model topology using a specified value of beta.
This script must be run before samples can be generated from the prior or
posterior.

#### cbayes_enumerate_sample.py ####

This script generates samples from the prior, or posterior, over model
topologies using the output a specified model probabilities file.

### GA Scripts ###


