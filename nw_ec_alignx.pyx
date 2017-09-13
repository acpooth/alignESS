# -*- coding: utf-8 -*-
#
#------------------------------
# Name:     nw_ec_align.px
# Purpose:  Needelman-Wunsh algorithm for EC number alignment -
#            Cython version
#
# @uthor:      acph - dragopoot@gmail.com
#
# Created:     Sep, 2017
# Copyright:   (c) acph 2017
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#------------------------------
"""Needelman-Wunsh algorithm for EC number alignment - Cython version.
This version performs ~30 times fastter than the original pure python/numpy
implementation."""

import cython
cimport cython
from cpython cimport bool
import numpy as np
cimport numpy as np
from libc.stdlib cimport rand
from libc.math cimport copysignf
cdef extern from "limits.h":
    int INT_MAX

DTYPE = np.int
ctypedef np.int_t DTYPE_t
DTYPEF = np.float
ctypedef np.float_t DTYPEF_t


cdef minarg(float f[3]):
    """Returns the minimum value and the index of that value from an array
    of 3 elements
    """
    cdef:
        int i
        int argmin = 0
        float mini = f[0]
    for i in range(1, 3):
        if f[i] < mini:
            mini = f[i]
            argmin = i
    return mini, argmin


# best performance
@cython.boundscheck(False)
@cython.wraparound(False)
def _FastNW(np.ndarray[DTYPEF_t, ndim=2] mat, dict ecs, list seq1, list seq2,
            float gap=0.9):
    """Perform dynamic programing alignment of EC numbers sequences. Creates
    the scoring and arrow matrices using brute force loops.

    Returns the dynamic programing matrix and the arrow matrix used to trace
    back the alignment.

    Arguments:
    - `mat`: EC number similarity (distance) matrix
    - `ecs`: EC number dictionary. Keys are EC numbers (3 levels of 
        classification) and values are indices to similarity matrix (mat)
    - `seq1`: enzymatic step sequence 1, list of (unaligned) EC numbers
    - `seq2`: enzymatic step sequence 2, list of (unaligned) EC numbers
    - `gap`: gap penalties, default = 1
    """
    cdef:
        int l1, l2, li, i, j, k, eci, ecj, minpos
        float mini, randum, randfactor
        float f[3]
    mini = 0
    l1, l2 = len(seq1), len(seq2)
    # Create the score and arrow matrices
    cdef np.ndarray[DTYPEF_t, ndim = 2] scoremat = np.zeros((l1 + 1, l2 + 1), DTYPEF)
    cdef np.ndarray[DTYPE_t, ndim = 2] arrow = np.zeros((l1 + 1, l2 + 1), DTYPE)
    # Create first row and first column with gaps
    for i in range(l2 + 1):
        scoremat[0, i] = i * gap
        arrow[0, i] = 1
    for i in range(l1 + 1):
        scoremat[i, 0] = i * gap
    # Fill the matrix. f array is cell posible values (len = 3)
    for i in range(1, l1 + 1):
        for j in range(1, l2 + 1):
            eci = ecs[seq1[i - 1]]
            ecj = ecs[seq2[j - 1]]
            sco = mat[eci, ecj]
            f[0] = scoremat[i - 1, j] + gap
            f[1] = scoremat[i, j - 1] + gap
            f[2] = scoremat[i - 1, j - 1] + sco
            # Adding some random variation to solve draws
            for k in range(3):
                randfactor = copysignf(0.001, f[k])
                randnum = rand() / float(INT_MAX)
                f[k] -= randfactor * randnum
            # select min value
            mini, minpos = minarg(f)
            # fill the matrices
            scoremat[i, j] = mini
            arrow[i, j] = minpos
    return arrow, mini


cdef gappen(list seq1, list seq2):
    """Returns the number
    Keyword Arguments:
    list seq1 -- 
    list seq2 -- 
    """
    cdef:
        int l1 = len(seq1)
        float gappen = 0
        str gap = '-.-.-'
        int i
        float sgaps = 0
        float ngaps = 0
        float sgaps2 = 0
        float ngaps2 = 0
        float fgaps = 0
        float fgaps2 = 0
        float fblock1 = 0
        float fblock2 = 0
        str ec1, ec2
        bool begin1 = True
        bool begin2 = True
        bool block1 = True
        bool block2 = True

    # four cases depending in wich sequence has gaps
    if gap not in seq1 and gap not in seq2:
        return gappen           # no gaps, return 0
    elif gap in seq1 and gap not in seq2:
        for i in range(l1):      # Gaps only in seq1
            ec1 = seq1[i]
            if begin1 == True:
                if ec1 == gap:
                    continue
                else:
                    begin1 = False
                    continue
            else:
                if ec1 == gap:
                    sgaps += 1
                    fblock1 += 1
                    if block1 == True:
                        ngaps += 1
                        block1 = False
                else:
                    block1 = True
                    fblock1 = 0
        if block1 == False:
            sgaps -= fblock1
            ngaps -= 1
        if sgaps == 0:
            gappen = 0.0
        else:
            gappen = ngaps / sgaps / 2
        return gappen
    elif gap not in seq1 and gap in seq2:
        for i in range(l1):      # Gaps only in seq2
            ec2 = seq2[i]
            if begin2 == True:
                if ec2 == gap:
                    continue
                else:
                    begin2 = False
                    continue
            else:
                if ec2 == gap:
                    sgaps += 1
                    fblock2 += 1
                    if block2 == True:
                        ngaps += 1
                        block2 = False
                else:
                    block2 = True
                    fblock2 = 0
        if block2 == False:
            sgaps -= fblock2
            ngaps -= 1
        if sgaps == 0:
            gappen = 0.0
        else:
            gappen = ngaps / sgaps / 2
    else:
        for i in range(l1):
            ec1 = seq1[i]
            ec2 = seq2[i]
            if begin1 == True:
                if ec1 == gap:
                    pass
                else:
                    begin1 = False
                    pass
            else:
                if ec1 == gap:
                    sgaps += 1
                    fblock1 += 1
                    if block1 == True:
                        ngaps += 1
                        block1 = False
                else:
                    block1 = True
                    fblock1 = 0
            ec2 = seq2[i]
            if begin2 == True:
                if ec2 == gap:
                    pass
                else:
                    begin2 = False
                    pass
            else:
                if ec2 == gap:
                    sgaps2 += 1
                    fblock2 += 1
                    if block2 == True:
                        ngaps2 += 1
                        block2 = False
                else:
                    block2 = True
                    fblock2 = 0
        if block1 == False:
            sgaps -= fblock1
            ngaps -= 1
        if block2 == False:
            sgaps2 -= fblock2
            ngaps2 -= 1
        if ngaps <= 0 and ngaps2 > 0:
            gappen = (ngaps2 / sgaps2) / 2
        elif ngaps > 0 and ngaps2 <= 0:
            gappen = (ngaps / sgaps) / 2
        elif ngaps <= 0 and ngaps2 <= 0:
            gappen = 0.0
        else:
            gappen = ((ngaps / sgaps) + (ngaps2 / sgaps2)) / 2
    return gappen

# better performance


@cython.boundscheck(False)
@cython.wraparound(False)
def scoring(np.ndarray[DTYPEF_t, ndim=2] mat, dict ecs, list seq1, list seq2,
            float gap=1, float fhomo=0.95, float fpengap=0.05):
    """ This function evaluates a pair of alignmed ESS and returns an entroy 
    based similarity (distance) score. This value is in the range 0-1. 0 means 
    more similar sequences and 1 means dissimilar sequences.

    The score tries to recall the fitness function proposed in Ortegon, et al. 
    2015. Comput Struct Biotechnol J. 9;13.

    Arguments:
    - `mat`: EC number similarity (distance) matrix
    - `ecs`: EC number dictionary. Keys are EC numbers (3 levels of 
        classification) and values are indices to similarity matrix (mat)
    - `seq1`: enzymatic step sequence 1, list of (unaligned) EC numbers
    - `seq2`: enzymatic step sequence 2, list of (unaligned) EC numbers
    - `gap`: gap penalties, default = 1
    - `fhomo`: Weight factor to multiplicate the homogeneity fraction of score
    - `fpengap`: Weight factor to multiplicate the gap penalization
    """
    cdef:
        unsigned int length, i, x, y
        float pengap, score
        float homo = 0.0
        str ec1, ec2
        str gapstr = '-.-.-'
        str dotstr = '...'
    # are sequences of same length?
    assert 0.99999 < (fhomo + fpengap) < 1.0001, \
        'fhomo + fpengap must sum aprox. 1'
    assert len(seq1) == len(seq2), \
        'The aligned sequences must be of the same length'
    length = len(seq1)  # Alignment length
    # -- Sum homogeneity evaluation
    for i in range(length):
        ec1 = seq1[i]
        ec2 = seq2[i]
        if ec1 == gapstr or ec2 == gapstr:  # and s1[i] != s2[1]
            homo += gap                     # gap homogeneity = 1
        elif ec1 == dotstr or ec2 == dotstr:
            pass
        elif ec1 == '' or ec2 == '':
            pass
        else:
            x = ecs[ec1]
            y = ecs[ec2]
            homo += mat[x, y]
    homo = homo / length        # mean homogeneity
    pengap = gappen(seq1, seq2)  # dap penalization
    # Final score calculation
    score = (homo * fhomo) + (pengap * fpengap)
    # there is no penalizarion for increment in columns
    return score


def _backtrace(np.ndarray[DTYPE_t, ndim=2] arrow, list seq1, list seq2):
    """Reads the arrow matrix of NW alignment and return two lists of
    with the aligned EC numbers.

    Arguments:
    - `arrow`: arrow matrix generated by dynamic programing
    - `seq1`: enzymatic step sequence, list of (unaligned) EC numbers
    - `seq2`: enzymatic step sequence, list of (unaligned) EC numbers
    """
    cdef:
        int v = len(seq1)
        int h = len(seq2)
        bool ok = True
        list st1 = []
        list st2 = []
    while ok:                 # backtrace walk
        if arrow[v, h] == 0:  # vertical best result, s1
            st1.append(seq1[v - 1])
            st2.append('-.-.-')
            v -= 1
        elif arrow[v, h] == 1:  # horizontal best result, s2
            st1.append('-.-.-')
            st2.append(seq2[h - 1])
            h -= 1
        elif arrow[v, h] == 2:  # diagonal best result, s1,s2 aligned
            st1.append(seq1[v - 1])
            st2.append(seq2[h - 1])
            v -= 1
            h -= 1
        if v == 0 and h == 0:
            ok = False
    # reverse sequences
    st1 = st1[::-1]
    st2 = st2[::-1]
    return st1, st2


def NW(np.ndarray[DTYPEF_t, ndim=2] mat, dict ecs, list seq1, list seq2,
       float gap=0.9, bool localize=False, bool nws=False,
       bool strfmt=True):
    """
    NW alignment function. Creates a pairwise alignment of Enzymatic Step 
    Sequences (ESS) using a Needelman-Wunsh algorithm.

    Arguments:
    - `mat`: EC number (3 levels) substitution matrix
    - `ecs`: List of ec numbers that represent the labels of the matrix
    - `seq1`: EC numbers sequence 1, list
    - `seq2`: EC numbers sequence 2, list
    - `gap`: gap penalty for NW algorithm (not for scoring)
    - `localize`: if localize = True, then, the function returns only the
            fragment of the alignment covered by the shortest sequence
            and the score is calculated accordingly
    - `nws`: if true, the function returns the scores as is returned by the
            NW algorithm
    - `strfmt`: if True. The aligned sequences are in string format. If 
            False, the sequences are in list format
    """
    cdef:
        float score, mini, scoring_gap
        str aseq1, aseq2
        int length, i, ii, fi
        char binit      # localize initial index boolena
        char bblock     # block mar boolean
    # create NW matrices
    arr, mini = _FastNW(mat, ecs, seq1, seq2, gap=gap)
    seq1, seq2 = _backtrace(arr, seq1, seq2)  # backtrace alingment
    # localized alignmet
    if localize:
        binit = 1
        bblock = 0
        length = len(seq1)
        ii = 0
        fi = 0
        for i in range(length):
            if seq1[i] != '-.-.-' and seq2[i] != '-.-.-' and binit:
                ii = i
                binit = 0
            elif (seq1[i] == '-.-.-' or seq2[i] == '-.-.-') and bblock:
                fi = i
                bblock = 0
            elif seq1[i] != '-.-.-' and seq2[i] != '-.-.-':
                bblock = 1
        if bblock:
            fi = length
        seq1 = seq1[ii: fi]
        seq2 = seq2[ii: fi]
    # NW score or entropy based score
    if nws:
        score = mini
    else:
        scoring_gap = 1
        score = scoring(mat, ecs, seq1, seq2, gap=scoring_gap, fhomo=0.95,
                        fpengap=0.05)
    if strfmt:
        aseq1 = ':'.join(seq1)
        aseq2 = ':'.join(seq2)
        return aseq1, aseq2, score
    else:
        return seq1, seq2, score
