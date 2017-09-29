#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:        nw_ec_align.px
# Purpose:     Program to align Enzimatic Step Sequences (ESS) using dynamic
#              programming
#
# @uthor:      acph - dragopoot@gmail.com
#
# Created:     Sep, 2017
# Copyright:   (c) acph 2017
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------
"""Program to align Enzimatic Step Sequences (ESS) using dynamic programming"""

import os
import argparse
import numpy as np
import pyximport
pyximport.install()
import nw_ec_alignx as nwx


# load similarity matrix
matfile = 'h_ent_mat.npz'       # ecs_entropy.py generated file
# exedir = os.path.dirname(argv[0])
# exepath = os.path.join(exedir, matfile)
exedir = os.path.realpath(__file__)
exedir = os.path.dirname(exedir)
exepath = os.path.join(exedir, matfile)

hmat = np.load(exepath)
ecs = hmat['ecs']
hmat = hmat['matrix']
lecs = list(ecs)
decs = {}
for i, v in enumerate(lecs):
    decs[v] = i

# -- End load matrices


def alignESS(seq1, seq2, gap):
    """Wrapper for nwx.NW function for ESS alignment

    Keyword Arguments:
    seq1 -- ESS 1. str, colon separated EC numbers.
    seq2 -- ESS 1. str, colon separated EC numbers.
    gap  -- Gap penalization.
    """
    seq1 = seq1.split(':')
    seq2 = seq2.split(':')
    seq1, seq2, score = nwx.NW(hmat, decs, seq1, seq2, gap=gap)
    return seq1, seq2, score


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ess1', type=str,
                        help='''ESS (3 levels EC numbers). Colon separated.
                        (1.2.3:3.5.-:...:9.9.9)''')
    parser.add_argument('ess2', type=str,
                        help='''ESS (3 levels EC numbers). Colon separated.
                        (1.2.3:3.5.-:...:9.9.9)''')
    parser.add_argument('--gap', type=float, default=0.9,
                        help='Gap penalization (from 0 to 1). Default = 0.9')
    args = parser.parse_args()
    return args


def main():
    from sys import stdout

    args = arg_parser()
    seq1, seq2, score = alignESS(args.ess1, args.ess2, gap=args.gap)

    stdout.write('ess1:\t{}\n'.format(seq1))
    stdout.write('ess2:\t{}\n'.format(seq2))
    stdout.write('score = {}\n'.format(score))


if __name__ == '__main__':
    main()
