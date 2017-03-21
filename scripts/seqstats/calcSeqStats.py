#!/usr/local/bin/python3
# calcSeqStats.py
# Calculate sequence statistics from FASTQ file.
# SUMMARY STATS: Number of sequences, read length, quality, G+C.
# DENSITY FILES: Used for density plots, etc.
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 23 Nov 2016
# Updated on 20 Mar 2017

from __future__ import absolute_import, division, print_function
import argparse
import os
import os.path
import sys
import sequtil as seq
import numpy as np
import util
import gzip


###############################################################################
## ARGUMENT PARSING
###############################################################################
parser = argparse.ArgumentParser(description="Calculate statistics from a "
                                 "FASTQ file.")
parser.add_argument("fastq", help="FASTQ file")
parser.add_argument("outdir", help="Directory for output files")
parser.add_argument("-s", "--summary_only", action="store_true",
                    help="Only produce a summary file")
parser.add_argument("--header", action="store_true",
                    help="Include header in output files")
parser.add_argument("--fasta", action="store_true",
                    help="File is in FASTA format")
parser.add_argument("--gzip", action="store_true",
                    help="File is compressed with GZIP")
parser.add_argument("--gc", action="store_true",
                    help="Calculate GC")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="Verbose output")

args = parser.parse_args()

seqFile = args.fastq
outdir = args.outdir if args.outdir.endswith("/") else args.outdir + "/"

# Check file exists
if not os.path.isfile(seqFile):
    util.printStatus("Sequence input file '" + fn + "' does not exist.")
    util.exitScript()

# Check output directory exists
if not os.path.isdir(outdir):
    util.printStatus("Output directory '" + outdir + "' does not exist.")
    util.exitScript()

# Give warning on file extension
sf_lc = seqFile.lower()
if args.gzip and not sf_lc.endswith(".gz"):
    util.printStatus("WARNING: file does not end in '.gz' - may not "
                     "be a GZIP file")
elif args.fasta and (not sf_lc.endswith(".fasta") and
                     not sf_lc.endswith(".fna") and
                     not sf_lc.endswith(".fa")):
    util.printStatus("WARNING: file does not end in '.fasta', '.fna',"
                     "nor '.fa' - may not be a FASTA file")
elif not sf_lc.endswith(".fastq") and sf_lc.endswith(".fq"):
    util.printStatus("WARNING: file does not end in '.fastq' not '.fq'"
                     "- may not be a FASTQ file")

fname = os.path.splitext(os.path.basename(seqFile))[0]
if args.gzip:
    # File name had two suffixes, must remove one more
    fname = os.path.splitext(fname)[0]

###############################################################################
## BEGIN PROCESSING FASTQ FILE
###############################################################################
# Intialize all data structures
numSequences = 0
allReadLengths = []
readLengths = {}
# Read positions will be used for dictionaries
readPositions = ["0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69",
                 "70-79", "80-89", "90-99", "100-109", "110-119", "120-129",
                 "130-139", "140-149", "150-159", "160-169", "170-179",
                 "180-189", "190-199", "200-209", "210-219", "220-229",
                 "230-239", "240-249", "250+"]
if not args.fasta:
    allqualities = {}
    qualities = {}

if args.gc:
    allgcs = {}
    gcs = {}


# Parse the file
# First determine if the file is FASTA or FASTQ format
if args.fasta:
    # Check if file is gzipped
    if args.gzip:
        fh = gzip.open(seqFile, "rt")
    else:
        fh = open(seqFile, "r")

    inSeq = False
    for l in fh:
        l = l.rstrip("\n")

        if l.startswith(">"):
            inSeq = True

        elif inSeq:
            # In sequence
            numSequences += 1
            if args.verbose:
                util.printStatus("Working on sequence " + str(numSequences),
                                 end="\r")
            # Add length of sequence read
            lseq = len(l)
            allReadLengths.append(lseq)
            if lseq not in readLengths:
                readLengths[lseq] = 0
            readLengths[lseq] += 1

            if args.gc:
                # Get current sequence GC percentage
                # Add to list of all GC percentages
                strgc = "{:.2f}".format(seq.GC(l))
                if strgc not in allgcs:
                    allgcs[strgc] = 0
                allgcs[strgc] += 1

                # Calculate average GC percentages per sequence interval
                # 10 bp intverals
                currGC = seq.GC_interval(l, 10)

                # Iterate through intervals and record GC percentages
                for idx, g in enumerate(currGC):
                    strgc = "{:.2f}".format(g)
                    # Get read interval String based on index
                    if idx >= len(readPositions):
                        rp = "250+"
                    else:
                        rp = readPositions[idx]

                    # Use GC percentage as key in dictionary
                    if rp not in gcs:
                        gcs[rp] = {}
                    if strgc not in gcs[rp]:
                        gcs[rp][strgc] = 0

                    gcs[rp][strgc] += 1

            # End of sequence read
            inSeq = False

else:
    # Check if file is gzipped
    if args.gzip:
        fh = gzip.open(seqFile, "rt")
    else:
        fh = open(seqFile, "r")

    inSeq = False
    inQual = False
    for l in fh:
        l = l.strip("\n")

        if l.startswith("@"):
            inSeq = True

        elif inSeq:
            # In sequence
            numSequences += 1
            if args.verbose:
                util.printStatus("Working on sequence " + str(numSequences),
                                 end="\r")
            # Add length of sequence read
            lseq = len(l)
            allReadLengths.append(lseq)
            if lseq not in readLengths:
                readLengths[lseq] = 0
            readLengths[lseq] += 1

            if args.gc:
                # Get current sequence GC percentage
                # Add to list of all GC percentages
                strgc = "{:.2f}".format(seq.GC(l))
                if strgc not in allgcs:
                    allgcs[strgc] = 0
                allgcs[strgc] += 1

                # Calculate average GC percentages per sequence interval
                # 10 bp intverals
                currGC = seq.GC_interval(l, 10)

                # Iterate through intervals and record GC percentages
                for idx, g in enumerate(currGC):
                    strgc = "{:.2f}".format(g)
                    # Get read interval String based on index
                    if idx >= len(readPositions):
                        rp = "250+"
                    else:
                        rp = readPositions[idx]

                    # Use GC percentage as key in dictionary
                    if rp not in gcs:
                        gcs[rp] = {}
                    if strgc not in gcs[rp]:
                        gcs[rp][strgc] = 0

                    gcs[rp][strgc] += 1

            # End of sequence read
            inSeq = False

        elif l.startswith("+"):
            # Next line will by quality scores
            inQual = True

        elif inQual:
            # Get current sequence quality scores
            # Add to list of all quality scores
            for q in seq.QualToInt(l):
                strq = "{:.2f}".format(q)
                if strq not in allqualities:
                    allqualities[strq] = 0
                allqualities[strq] += 1

            # Calculate average quality values per sequence inverval
            # 10bp intervals
            currQ = seq.QualToInt_interval(l)

            # Iterate through intervals and record quality values
            for idx, q in enumerate(currQ):
                strq = "{:.2f}".format(q)
                # Get read interval String based on index
                if idx >= len(readPositions):
                    rp = "250+"
                else:
                    rp = readPositions[idx]

                # Use quality score as key in dictionary
                if rp not in qualities:
                    qualities[rp] = {}
                if strq not in qualities[rp]:
                    # Initialize counter if does not exist
                    qualities[rp][strq] = 0

                qualities[rp][strq] += 1

            # End of quality scores
            inQual = False
fh.close()
# End file parsing
print("", file=sys.stderr)

###############################################################################
## CALCULATE AVERAGES
###############################################################################
# Values for summary file
meanLen = np.mean(allReadLengths)
stdevLen = np.std(allReadLengths)
if not args.fasta:
    quals = list(allqualities.keys())
    counts = [allqualities[x] for x in quals]
    meanQual = np.average(float(quals), weights=counts)
    stdevQual = np.sqrt(np.average((float(quals) - meanQual) ** 2,
                                   weight=counts))
if args.gc:
    gcpercs = list(allgcs.keys())
    counts = [allgcs[x] for x in gcpercs]
    meanGC = np.average(float(gcpercs), weights=counts)
    stdevGC = np.sqrt(np.average((float(gcpercs) - meanGC) ** 2,
                                 weight=counts))

# Quality score averages
if not args.fasta:
    qavgs = {}  # Store read positions (key) and their average (value)
    for rp, qscores in qualities.items():
        n = 0  # Number of instances in that read position
        total = 0  # Total Q score
        for q, count in qscores.items():
            n += count
            total += float(q) * count
        # Calculate average quality score for that read postion
        qavgs[rp] = total / n

# GC percentage averages
if args.gc:
    gcavgs = {}  # Store read positions (key) and their average (value)
    for rp, gcpercs in gcs.items():
        n = 0  # Number of instances in that read position
        total = 0  # Total GC percentage
        for gc, count in gcpercs.items():
            n += count
            total += float(gc) * count
        # Calculate average GC percentage for that read position
        gcavgs[rp] = total / n

###############################################################################
## BEGIN OUTPUT FILES
###############################################################################
if args.verbose:
    util.printStatus("Printing summary file")
with open(outdir + fname + "_summary.tsv", "w") as f:
    if args.header:
        f.write("Metric\tValue\n")
    f.write("Number of sequences\t{}\n".format(numSequences))
    f.write("Read length (mean)\t{:.2f}\n".format(meanLen))
    f.write("Read length (standard deviation)\t{:.2f}\n".format(stdevLen))
    if not args.fasta:
        f.write("Quality (mean)\t{:.2f}\n".format(meanQual))
        f.write("Quality (standard deviation)\t{:.2f}\n".format(stdevQual))
        for rp in readPositions:
            if rp not in qualities:
                break
            f.write("Quality by interval ({})\t".format(rp))
            f.write("{:.2f}\n".format(qavgs[rp]))
    if args.gc:
        f.write("GC ratio (mean)\t{:.2f}\n".format(meanGC))
        f.write("GC ratio (standard deviation)\t{:.2f}\n".format(stdevGC))
        for rp in readPositions:
            if rp not in gcs:
                break
            f.write("GC ratio by interval ({})\t".format(rp))
            f.write("{:.2f}\n".format(gcavgs[rp]))

# Print out density files if necessary
if not args.summary_only:
    if args.verbose:
        util.printStatus("Printing read lengths file")
    with open(outdir + fname + "_readlengths.tsv", "w") as f:
        if args.header:
            f.write("Length\tCount\n")
        for l, count in readLengths.items():
            f.write("{}\t{}\n".format(l, count))

    if not args.fasta:
        if args.verbose:
            util.printStatus("Printing qualities file")
        with open(outdir + fname + "_qualities.tsv", "w") as f:
            if args.header:
                f.write("Position\tQuality\tCount\n")
            for rp in readPositions:
                if rp not in qualities:
                    break
                for q, count in qualities[rp].items():
                    f.write("{}\t{}\t{}\n".format(rp, q, count))
        with open(outdir + fname + "_allqualities.tsv", "w") as f:
            if args.header:
                f.write("Quality\tCount\n")
            for q, count in allqualities.items():
                f.write("{:.2f}\t{}\n".format(q, count))

    if args.gc:
        if args.verbose:
            util.printStatus("Printing GC file")
        with open(outdir + fname + "_gcratios.tsv", "w") as f:
            if args.header:
                f.write("Position\tGC_ratio\tCount\n")
            for rp in readPositions:
                if rp not in gcs:
                    break
                for gc, count in gcs[rp].items():
                    f.write("{}\t{}\t{}\n".format(rp, gc, count))
        with open(outdir + fname + "_allgcratios.tsv", "w") as f:
            if args.header:
                f.write("GC_ratio\tCount\n")
            for gc, count in allgcs.items():
                f.write("{:.2f}\t{}\n".format(gc, count))
