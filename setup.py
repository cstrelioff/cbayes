#!/usr/bin/env python

from distutils.core import setup

setup(name='cbayes',
    version='0.1',
    description='Scripts for Bayesian inference of uHMMs and eMs using CMPy.',
    author='Christopher C. Strelioff',
    author_email='chris.strelioff@gmail.com',
    url='https://github.com/cstrelioff/cbayes',
    packages = ['cbayes'],
    scripts = ['bin/cbayes_bash_enumerate_convergence.py',
        'bin/cbayes_bash_enumerate_overlap_analysis.py',
        'bin/cbayes_gabayes.py',
        'bin/cbayes_enumerate_AddToDB.py',
        'bin/cbayes_enumerate_CalcProbs.py',
        'bin/cbayes_enumerate_PriorAddToDB.py',
        'bin/cbayes_enumerate_ProcessSamples.py',
        'bin/cbayes_enumerate_Sample.py',
        'bin/cbayes_slurm_enumerate_convergence.py',
        'bin/cbayes_slurm_enumerate_overlap_analysis.py',
        'bin/cbayes_slurm_gabayes.py',
        'bin/cbayes_create_process_datafile.py']
    )
