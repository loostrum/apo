#!/bin/bash
# Script for processing eShel spectra
# Depends on specnorm.py and order_merging.py
#
# Author: Leon Oostrum
# E-Mail: l.c.oostrum@uva.nl

trap exit INT

self=`readlink -f $0`
scriptdir=`dirname $self`
specnorm=$scriptdir/specnorm.py
merging=$scriptdir/order_merging.py

# Check if specnorm is available
if [[ ! -f $specnorm ]]; then
    echo "specnorm.py not found."
    echo "Please make sure it is located in the same directory as this script."
    exit 1
fi

# Check if order merging is available
if [[ ! -f $merging ]]; then
    echo "order_merging.py not found."
    echo "Please make sure it is located in the same directory as this script."
    exit 1
fi

# Check if data directory is supplied by user
if [[ "$#" -eq "0" ]]; then
    echo "Usage: `basename $0` DATADIR SAVEDIR"
    echo "DATADIR: Directory containing the fully calibrated files (1C) produced by Audela, one per order."
    echo "SAVEDIR: Optional, default is DATADIR/../norm. Normalized spectra are saved to this directory."
    exit 1
fi

# set directories
datadir=`readlink -f $1`
if [[ ! -d $datadir ]]; then
    echo "$datadir does not exist."
    exit 1
fi

# Check if SAVEDIR is specified
# Set to default if not specified, then ask if this is ok.
if [[ "$#" -eq "1" ]]; then
    # Set default and ask if this is ok
    savedir=`readlink -f $datadir/../norm`
    echo "SAVEDIR not specified. Set to $savedir"
    read -rep "Is this ok [Y/n]?" ans
    case $ans in 
        [Nn] ) exit 1;;
        [Yy] ) ;;
        *    ) ;;
    esac
else
    savedir=`readlink -f $2`
fi

# Check if savedir exists
if [[ ! -d "$savedir" ]]; then
    echo "$savedir does not exist."
    exit 1
fi

# Create list of files
files=`ls $datadir/*_1C_[0-9][0-9].fit 2>/dev/null`
# There has to be better way to do this ...
nfiles=`ls $datadir/*_1C_[0-9][0-9].fit 2>/dev/null | wc -l`

# check if directory contains the right files
if [[ "$nfiles" -eq "0" ]]; then
    echo "Could not find files."
    echo "Did you select the right directory?"
    exit 1
fi


# Call the actual normalization script for each file
echo "Found $nfiles files."
echo "Will now call normalization script."

for file in $files; do
    echo `basename $file`
    $specnorm $file $savedir 
done

# Merge orders
echo "Normalization finished."
read -rep "Start order merging [Y/n]?" ans
case $ans in 
    [Nn] ) exit 0;;
    [Yy] ) $merging $savedir $savedir;;
    *    ) $merging $savedir $savedir;;
esac

# We are done :D
exit 0
