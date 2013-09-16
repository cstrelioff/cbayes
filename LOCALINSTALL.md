# Local Install Notes #

The scripts in `cmpy_bayes` are intended to be used as system scripts. In order
for this to work the scripts have to be in the system path. In general this will
be an *OS-specific* issue, so we will separate these out by OS.

## Ubuntu ##

Installation with the `--user` switch will install to `.local/bin` in the user
home directory on Ubuntu.  This location is **not** in the system path by
default, so we add this location to the path.  In the `.bashrc` file (create a
file in your home directory if it does not already exist) add the following
lines:

    # include .local/bin for local python scripts
    export PATH=~/.local/bin/:$PATH

To make this definition active for the current session you have to `source` the
`.bashrc` file using:

    source .bashrc

Turns out that an interactive login does not source `.bashrc`, for example when
you want to use these scripts via `ssh` on a remote server.  So, create a
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

## MAC ##

No notes.

## Windows ##

No notes.
