#!/bin/bash
# seqstats.sh
# Sequence statistics pipeline
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 23 Nov 2016
# Updated on 23 Nov 2016

VERSION="0.1"

usage() {
    echo "$1
$scriptname version $VERSION

usage: $scriptname -f fastq -d output_dir [Options]

Required
   -f [fastq]              : FASTQ file
   -o [output_dir]         : Directory for output files

Optional
   -t                      : Title for plots
   -v                      : Verbose output
   -h, -?, --help          : This help message

Notes
   - None

" >&2
}


error() {
    echo "*****FATAL ERROR OCCURRED*****" >&1 | tee -a $log
    echo $1 >&1 | tee -a $log
    exit 1
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
fastq=""
outdir=""
title=""
verbose=0

# Set pipefail for catching errors in piped commands
set -o pipefail

while [[ $# != 0 ]]; do
    case $1 in
    -h|-\?|--help)
        usage
        exit 2
        ;;
    -f)
        shift
        [[ ! $1 || $(printf "%s" "$1" | perl -ne 'm/(^-.$)/; print $1;') ]] && echo "Missing -f value" >&2 && usage && exit 2
        fastq=$1
        ;;
    -o)
        shift
        [[ ! $1 || $(printf "%s" "$1" | perl -ne 'm/(^-.$)/; print $1;') ]] && echo "Missing -o value" >&2 && usage && exit 2
        outdir=$1
        ;;
    -t)
        shift
        [[ ! $1 || $(printf "%s" "$1" | perl -ne 'm/(^-.$)/; print $1;') ]] && echo "Missing -t value" >&2 && usage && exit 2
        title=$1
        ;;
    -v)
        verbose=1
        ;;
    *)
        echo "Unknown option $1" >&2
        usage
        exit 2
    esac
    shift
done

# Check if required variables are set
if [[ ! $fastq || ! $outdir ]]; then
    usage "Missing one or more required arguments."
    exit 2
fi

# Extract name from FASTQ file
name=$(basename $fastq)
name=${name%.*}

# Check if title was given
if [[ ! $title ]]; then
    title=$name
fi

# Begin statistics scripts
getTime && echo "${currtime}    *****Starting sequence statistics scripts*****"  >&1
if (( !$verbose )); then
    getTime && echo "${currtime}    Note: verbose flag was not set."  >&1
fi

# Calculate sequencing stats
cmd="python2.7 calcSeqStats.py $fastq $outdir --header"
(( $verbose )) && getTime && echo "${currtime}    Executing $cmd"  >&1
eval $cmd  2>&1
[[ $? -ne 0 ]] && getTime && error "${currtime}    Fail on command: $cmd"

# Plot quality stats
cmd="Rscript seqstats_density.R -i ${outdir}/${name}_qualities -d $outdir --header -s qualities -t $title"
(( $verbose )) && getTime && echo "${currtime}    Executing $cmd"  >&1
eval $cmd  2>&1
[[ $? -ne 0 ]] && getTime && error "${currtime}    Fail on command: $cmd"

# Plot GC ratio stats
cmd="Rscript seqstats_density.R -i ${outdir}/${name}_gcratios -d $outdir --header -s gcratios -t $title"
(( $verbose )) && getTime && echo "${currtime}    Executing $cmd"  >&1
eval $cmd  2>&1
[[ $? -ne 0 ]] && getTime && error "${currtime}    Fail on command: $cmd"

# Plot sequence length stats
############################
## SKIPPING ##
#cmd="Rscript seqstats_readLengths.R -i ${outdir}/${name}_readlengths -d $outdir --header -t $title"
#(( $verbose )) && getTime && echo "${currtime}    Executing $cmd"  >&1
#eval $cmd  2>&1
#[[ $? -ne 0 ]] && getTime && error "${currtime}    Fail on command: $cmd"

getTime && echo "${currtime}    *****Completed!*****"  >&1
