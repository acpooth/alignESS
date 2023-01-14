#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:
# Purpose:
#
# @uthor:      acph - dragopoot@gmail.com
#
# Created:
# Copyright:   (c) acph 2016
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------
"""DON't use Old utility script for recreating alignments for a list of EC ids and socores"""

outf = open('alignments_localize.txt', 'w')

for i, j, s in res:
    ess1 = glyimp[0][i-1]
    j_ = nrseq[1].index(str(j))
    ess2 = nrseq[0][j_]
    ea1, ea2, s = alignESS.nwx.NW(hmat, decs,
                                  ess1, ess2,
                                  localize=True, strfmt=True)
    ali = (f'{i}\t{ea1}\n'
           f'{j}\t{ea2}\n'
           f'score: {s}\n///\n\n')
    outf.write(ali)

outf.close()


outf = open('alignments_complete.txt', 'w')

for i, j, s in res:
    ess1 = glyimp[0][i-1]
    j_ = nrseq[1].index(str(j))
    ess2 = nrseq[0][j_]
    ea1, ea2, s = alignESS.nwx.NW(hmat, decs,
                                  ess1, ess2,
                                  localize=False, strfmt=True)
    ali = (f'{i}\t{ea1}\n'
           f'{j}\t{ea2}\n'
           f'score: {s}\n///\n\n')
    outf.write(ali)

outf.close()
