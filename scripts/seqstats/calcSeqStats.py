#!/usr/local/bin/python2.7
# calcSeqStats.py
# Calculate sequence statistics from FASTQ file.
# SUMMARY STATS: Number of sequences, read length, quality, G+C.
# DENSITY FILES: Used for density plots, etc.
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 23 Nov 2016
# Updated on 28 Nov 2016

from __future__ import absolute_import, division, print_function
import argparse
import os
import os.path
import sys
import sequtil as seq
import numpy as np
import util


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

args = parser.parse_args()

fastq = args.fastq
outdir = args.outdir if args.outdir.endswith("/") else args.outdir + "/"

# Check file exists
if not os.path.isfile(fastq):
    util.printStatus("FASTQ file '" + fn + "' does not exist.")
    util.exitScript()

# Check output directory exists
if not os.path.isdir(outdir):
    util.printStatus("Output directory '" + outdir + "' does not exist.")
    util.exitScript()

# Give warning on file extension
if not fastq.lower().endswith(".fastq"):
    util.printStatus("WARNING: file does not end in '.fastq' - may not "
          "be a FASTQ file")

fname = os.path.splitext(os.path.basename(fastq))[0]

###############################################################################
## BEGIN PROCESSING FASTQ FILE
###############################################################################
# Intialize all data structures
numSequences = 0
readLengths = []
readPositions = ["0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69",
                 "70-79", "80-89", "90-99", "100-109", "110-119", "120-129",
                 "130-139", "140-149", "150-159", "160-169", "170-179",
                 "180-189", "190-199", "200-209", "210-219", "220-229",
                 "230-239", "240-249", "250+"]
allqualities = []
qualities = {k: [] for k in readPositions}
allgcs = []
gcs = {k: [] for k in readPositions}


# Parse FASTQ
with open(fastq, "r") as fh:
    inSeq = False
    inQual = False
    for l in fh:
        l = l.strip()

        if l.startswith("@"):
            inSeq = True

        elif inSeq:
            # In sequence
            numSequences += 1
            readLengths.append(len(l))
            allgcs.append(seq.GC(l))
            currGC = seq.GC_interval(l, 10)
            for idx, g in enumerate(currGC):
                if idx >= len(readPositions):
                    rp = "250+"
                else:
                    rp = readPositions[idx]
                gcs[rp].append(g)
            inSeq = False

        elif l.startswith("+"):
            inQual = True

        elif inQual:
            allqualities.append(seq.QualToInt(l))
            currQ = seq.QualToInt_interval(l)
            for idx, q in enumerate(currQ):
                if idx >= len(readPositions):
                    rp = "250+"
                else:
                    rp = readPositions[idx]
                qualities[rp].append(q)
            inQual = False
# End file parsing

###############################################################################
## BEGIN OUTPUT FILES
###############################################################################
# Values for summary file
meanLen = np.mean(readLengths)
stdevLen = np.std(readLengths)
meanQual = np.mean(allqualities)
stdevQual = np.std(allqualities)
meanGC = np.mean(allgcs)
stdevGC = np.std(allgcs)
with open(outdir + fname + "_summary.txt", "w") as f:
    if args.header:
        f.write("Metric\tValue\n")
    f.write("Number of sequences\t{}\n".format(numSequences))
    f.write("Read length (mean)\t{:.3f}\n".format(meanLen))
    f.write("Read length (standard deviation)\t{:.3f}\n".format(stdevLen))
    f.write("Quality (mean)\t{:.3f}\n".format(meanQual))
    f.write("Quality (standard deviation)\t{:.3f}\n".format(stdevQual))
    for rp in readPositions:
        if len(qualities[rp]) == 0:
            break
        f.write("Quality by interval ({})\t".format(rp))
        f.write("{:.3f}\n".format(np.mean(qualities[rp])))
    f.write("GC ratio (mean)\t{:.3f}\n".format(meanGC))
    f.write("GC ratio (standard deviation)\t{:.3f}\n".format(stdevGC))
    for rp in readPositions:
        if len(gcs[rp]) == 0:
            break
        f.write("GC ratio by interval ({})\t".format(rp))
        f.write("{:.3f}\n".format(np.mean(gcs[rp])))

# Print out density files if necessary
if not args.summary_only:
    with open(outdir + fname + "_readlengths", "w") as f:
        if args.header:
            f.write("Lengths\n")
        for l in readLengths:
            f.write(str(l) + "\n")

    with open(outdir + fname + "_qualities", "w") as f:
        if args.header:
            f.write("Position\tQuality\n")
        for rp in readPositions:
            for q in qualities[rp]:
                f.write("{}\t{:.3f}\n".format(rp, q))

    with open(outdir + fname + "_gcratios", "w") as f:
        if args.header:
            f.write("Position\tGC_ratio\n")
        for rp in readPositions:
            for gc in gcs[rp]:
                f.write("{}\t{:.3f}\n".format(rp, gc))
