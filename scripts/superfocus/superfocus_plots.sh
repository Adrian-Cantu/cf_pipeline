#!/bin/bash
# superfocus_plots.sh
# Create figures for SUPER-FOCUS output
# GitHub: https://github.com/metageni/SUPER-FOCUS/
# Ref: Silva GGZ, Green K., B. E. Dutilh, and R. A. Edwards: SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data.
#      (Bioinformatics. 2015 Oct 9. pii: btv584.)
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 08 Nov 2016
# Updated on 08 Nov 2016

VERSION="0.1"

usage() {
    echo "$1
usage: $scriptname -d SF_dir [Options]

Required
   -d [SF_dir]             : SUPER-FOCUS directory of files

Optional
   -h, -?, --help          : This help message
   -v                      : Verbose output

" >&2
}


getTime() {
    currtime=$(date "+[%F %H:%M:%S]")
}


timeStamp() {
    timestamp=$(date "+%Y%m%dT%H%M%S")
}


####################################################
#ARGUMENT PARSING
####################################################
scriptname=$(echo $0 | perl -ne '/\/?.*\/(.+)/; print $1;')
sfdir=""
verbose=""

# Set pipefail for catching errors in piped commands
set -o pipefail

while [[ $# != 0 ]]; do
    case $1 in
    -h|-\?|--help)
        usage $frmmax
        exit 2
        ;;
    -d)
        shift
        [[ ! $1 || $(printf "%s" "$1" | perl -ne 'm/(^-.$)/; print $1;') ]] && echo "Missing -d value" >&2 && usage && exit 2
        sfdir=$1
        ;;
    -v)
        verbose=0;
        ;;
    *)
        echo "Unknown option $1" >&2
        usage
        exit 2
    esac
    shift
done

# Check if required variables are set
if [[ ! $sfdir ]]; then
    usage "Missing one or more required arguments."
    exit 2
fi


# Create log file
timeStamp
log="log_${scriptname}_${timestamp}.txt"

# Plotting functions
getTime && echo "${currtime}    *****Starting plotting scripts!*****"  >&1 | tee -a $log
cmd="Rscript superfocus_functions.R -d ${sfdir}/"
[[ $verbose ]] && getTime && echo "${currtime}    Executing $cmd"  >&1 | tee -a $log
eval $cmd  2>&1 | tee -a $log
[[ $? -ne 0 ]] && getTime && error "${currtime}    Fail on command: $cmd"

# Plotting virulence functions
cmd="Rscript superfocus_virulence.R -d ${sfdir}/"
[[ $verbose ]] && getTime && echo "${currtime}    Executing $cmd"  >&1 | tee -a $log
eval $cmd  2>&1 | tee -a $log
[[ $? -ne 0 ]] && getTime && error "${currtime}    Fail on command: $cmd"

getTime && echo "${currtime}    *****Completed!*****"  >&1 | tee -a $log
