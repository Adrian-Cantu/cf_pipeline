# sequtil.py
# Utility functions for sequence data
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com
# Created on 23 Nov 2016
# Updated on 23 Nov 2016

from __future__ import absolute_import, division, print_function
import util
import numpy as np


def GC(seq):
    """Calculate GC ratio in given sequence"""
    # Check that sequence is non-empty
    seqlen = len(seq)
    if seqlen == 0:
        util.printStatus("WARNING in GCinterval(): sequence is an "
                         "empty string")
        return None

    seq = seq.lower()
    return (seq.count("g") + seq.count("c")) / seqlen


def GC_interval(seq, interval=10):
    """Calculate GC ratio in given sequence per interval"""
    # Check that interval is positive
    if interval < 1:
        util.printStatus("WARNING in GC_interval(): cannot use an interval "
                         "less than 1. Defaulting to interval = 10")
        interval = 10

    # Check that sequence is non-empty
    seqlen = len(seq)
    if seqlen == 0:
        util.printStatus("WARNING in GC_interval(): sequence is an "
                         "empty string")
        return None

    seq = seq.lower()
    gcs = []
    for i in xrange(0, seqlen, interval):
        currSeq = seq[i: i + interval]
        cslen = len(currSeq)
        gcs.append((currSeq.count("g") + currSeq.count("c")) / cslen)

    return gcs


def QualToInt(qual, phred=33):
    """Convert quality character to integer value"""
    # Check that quality is non-empty
    qlen = len(qual)
    if qlen == 0:
        util.printStatus("WARNING in QualToInt(): quality is an "
                         "empty string")
        return None

    # Check that quality is either 33 or 64
    if not phred == 33 and not phred == 64:
        util.printStatus("WARNING in QualToInt(): phred value is not the "
                         "expected 33 or 64")

    quals = []
    for nt in qual:
        quals.append(ord(nt) - phred)

    return quals


def QualToInt_interval(qual, phred=33, interval=10):
    """Convert quality character to integer value and calculate average"""
    # Check that interval is positive
    if interval < 1:
        util.printStatus("WARNING in QualToInt_interval(): cannot use an "
                         "interval less than 1. Defaulting to interval = 10")
        interval = 10

    # Check that quality is non-empty
    qlen = len(qual)
    if qlen == 0:
        util.printStatus("WARNING in QualToInt_interval(): quality is an "
                         "empty string")
        return None

    # Check that quality is either 33 or 64
    if not phred == 33 and not phred == 64:
        util.printStatus("WARNING in QualToInt_interval(): phred value is not "
                         "the expected 33 or 64")

    quals = []
    for i in xrange(0, qlen, interval):
        currQuals = qual[i: i + interval]
        quals.append(np.mean([ord(nt) - phred for nt in currQuals]))

    return quals
