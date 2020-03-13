#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
# Name:     test_scores-nw.py
# Purpose:  Test scores in old and new nw_ec_align
#
# @uthor:      acph - dragopoot@gmail.com
#
# Created:     Sep, 2017
# Copyright:   (c) acph 2017
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#------------------------------
"""Test scores in old and new nw_ec_align"""

import nw_ec_align as nw
lecs = nw.lecs
decs = nw.decs

import pyximport
pyximport.install()
import dlist
import nw_ec_alignx as nwx


s1 = '6.2.1:2.3.1:1.2.4:2.7.1:4.2.1:5.4.2:2.7.2:1.2.1:4.1.2:2.7.1:5.3.1:5.4.2:3.1.3'
s1_ = '2.3.1:1.2.4:2.7.1:5.4.2:2.7.2:1.2.1'
s2 = '1.2.4:2.3.1:2.3.3:4.1.1'
s3 = '6.3.3:6.3.4:5.4.99:6.3.2:4.3.2:2.4.2:3.2.2'
ss1 = s1.split(':')
ss1_ = s1_.split(':')
ss2 = s2.split(':')
ss3 = s3.split(':')

ess = [ss1, ss2, ss3, ss1_]
print("Test ESS:")
print('-'*45)
print('1\t', nw.s1)
print('2\t', nw.s2)
print('3\t', nw.s3)
print('4\t', nw.s1_)
print('-'*45, '\n')
print('\tSCORES')
print('-'*45)
print('\t|     Old\t|    New')
print('s1-s2\t|Sco\tLoc\t|Sco\tLoc')
print('-'*45)

line = '{}-{}\t|{:.3}\t{:.3}\t|{:.3}\t{:.3}'
for ii, i in enumerate(ess):
    for ji, j in enumerate(ess):
        sold = nw.NW(nw.hmat, lecs, ':'.join(i), ':'.join(j))[2]
        snew = nwx.NW(nw.hmat, decs, i, j)[2]
        soldl = nw.NW(nw.hmat, lecs, ':'.join(i), ':'.join(j),
                      local=True, localize=True)[2]
        snewl = nwx.NW(nw.hmat, decs, i, j, localize=True)[2]
        print(line.format(ii, ji, sold, soldl, snew, snewl))



