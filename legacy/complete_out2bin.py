#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
# Name:     
# Purpose:  
# 
# @uthor:   acph - dragopoot@gmail.com
#
# Created:     
# Copyright:   (c) acph 2014
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#------------------------------
"""
WARNING - Text arrays are to large, try DB
Transform the text output of alignAll_complete_out.py into a numpy
simple binary. The output file is tab delimited.

The binary contains 6 arrays.
ids1
ids2
scores
ess1
ess2

use:
complete_out2bin.py OUT_FILE

"""

from sys import argv
from subprocess import Popen, PIPE
import numpy as np


def create_arrays(in_file, sco_dtype=np.float32):
    """
    Returns the empty arrays with the corresponding length.
    Uses GNU wc.
    """
    print "Creating arrays"
    cline = "wc -l {}".format(in_file).split()
    child = Popen(cline, stdout=PIPE)
    lines = child.stdout.read()
    lines = int(lines.split()[0])
    ids1 = np.zeros(lines, dtype=np.uint16)
    ids2 = np.zeros(lines, dtype=np.uint16)
    scores = np.zeros(lines, dtype=sco_dtype)
    ess1 = np.zeros(lines, dtype='S210')
    ess2 = np.zeros(lines, dtype='S210')
    return ids1, ids2, scores, ess1, ess2

def fill_arrays(in_file, sco_dtype=np.float32):
    """Fills the arrays created by create_arrays for the file in_file
    and save in a numpy compressed file
    
    Arguments:
    - `in_file`: Complete output of the all vs all alignments
    - `sco_dtype`: set the datafile of the scores
    """
    ids1, ids2, scores, ess1, ess2 = create_arrays(in_file,
                                                   sco_dtype=sco_dtype)
    i = 0
    print "Filling arrays ..."
    with open(in_file) as inf:
        for line in inf:
            line = line.strip()
            i1, i2, sco, s1, s2 = line.split('\t')
            ids1[i] = i1
            ids2[i] = i2
            scores[i] = sco
            ess1[i] = s1
            ess2[i] = s2
            i += 1
    outfname = in_file + '.npz' 
    np.savez_compressed(outfname, ids1=ids1, ids2=ids2, scores=scores,
                        ess1=ess1, ess2=ess2)
    print "  :)  Done!!!"

def main():
    """
    """
    try:
        infile = argv[1]
    except:
        print "First (only) argument must be the complete output file"
        exit()
    fill_arrays(infile)

if __name__ == '__main__':
    main()