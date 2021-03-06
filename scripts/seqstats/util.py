# util.py
# Utility functions
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com
# Created on 23 Nov 2016
# Updated on 13 Mar 2017

from __future__ import absolute_import, division, print_function
import sys
import time
import datetime


def timeStamp():
    '''Return time stamp'''
    t = time.time()
    fmt = '[%Y-%m-%d %H:%M:%S]'
    return datetime.datetime.fromtimestamp(t).strftime(fmt)


def printStatus(msg, end="\n"):
    '''Print status message'''
    print('{}    {}'.format(timeStamp(), msg), file=sys.stderr, end=end)
    sys.stderr.flush()


def exitScript(num=1):
    '''Exit script'''
    sys.exit(num)
