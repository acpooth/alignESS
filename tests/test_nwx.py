#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     test_nwx.py
# Purpose:  Testing of nw_ec_alignx and alignESS modules
#
# @uthor:      acph - dragopoot@gmail.com
#
# Created:     Sep 2017, after shock!
# Copyright:   (c) acph 2017
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------
"""Testing of nw_ec_alignx and alignESS modules"""

import os
import sys
import pytest
from random import choice
# Fixing path for pytest
mod_path = os.path.abspath(__file__)
mod_path = os.path.dirname(mod_path)
sys.path.insert(0, mod_path + '/../')
# --End Fixing path for pytest
pyximport = pytest.importorskip('pyximport')
alignESS = pytest.importorskip('alignESS')
try:
    import nw_ec_align as nw
except ImportError:
    pass
pyximport.install()
import nw_ec_alignx as nwx

# use fixtures instead these commands
# decs = alignESS.decs
# hmat = alignESS.hmat


@pytest.fixture(scope="module")
def seqs_s():
    import sqlite3 as s3
    with s3.connect('tests/nr_part.db') as db:
        x = db.execute("SELECT ec3 from nrseqs")
        seqs_s = [i[0] for i in x]
    return seqs_s


@pytest.fixture(scope="module")
def seqs(seqs_s):
    seqs = [i.split(':') for i in seqs_s]
    return seqs


@pytest.fixture(scope="module")
def decs():
    return alignESS.decs


@pytest.fixture(scope="module")
def lecs():
    return alignESS.lecs


@pytest.fixture(scope="module")
def hmat():
    return alignESS.hmat


###########
# Markers #
###########
skipif = pytest.mark.skipif
xfail = pytest.mark.xfail
oldnw = pytest.mark.skipif('nw_ec_align' not in sys.modules,
                           reason='Requires nw_ec_align legacy version')

skipifnotfiles = skipif((not os.path.exists('tests/random_frag_ind.txt') or
                         not os.path.exists('tests/random_frag.txt') or
                         not os.path.exists('tests/nr_part.db') or
                         not os.path.exists('tests/part_name.txt')),
                        reason="""At least one test file is not present:
./tests/random_frag.txt
        random_frag_ind.txt
        nr_part.db
        part_name.txt """)

#####################
# Testing algorithm #
#####################


@oldnw
def test_scoring(seqs, seqs_s, decs, lecs, hmat):
    """Test for similarity in old and new alignmet scoring
    """

    assert len(seqs) == len(seqs_s)
    n = len(seqs)
    tests = 100
    for k in range(tests):
        i = choice(range(n))
        j = choice(range(n))
        s1, s2 = seqs[i], seqs[j]
        arr = nwx._FastNW(hmat, decs, s1, s2)
        # Cython aling and scoring
        as1, as2 = nwx._backtrace(arr[0], s1, s2)
        as1s = ':'.join(as1)
        as2s = ':'.join(as2)
        sco = nwx.scoring(hmat, decs, as1, as2)
        scos = nw.scoring(hmat, lecs, as1s, as2s)
        assert abs(sco - scos) <= 1e-6


@oldnw
def test_NW_score(hmat, decs, lecs):
    """Test for similarity in old and new alignment scoring
    Complete NW version
    """
    s1 = '6.2.1:2.3.1:1.2.4:2.7.1:4.2.1:5.4.2:2.7.2:1.2.1:4.1.2:2.7.1:5.3.1:5.4.2:3.1.3'
    s1_ = '2.3.1:1.2.4:2.7.1:5.4.2:2.7.2:1.2.1'
    # s2 = '1.2.4:2.3.1:2.3.3:4.1.1'
    # s3 = '6.3.3:6.3.4:5.4.99:6.3.2:4.3.2:2.4.2:3.2.2'
    ss1 = s1.split(':')
    ss1_ = s1_.split(':')
    # ss2 = s2.split(':')
    # ss3 = s3.split(':')
    ess = [ss1, ss1_]
    for ii, i in enumerate(ess):
        for ji, j in enumerate(ess):
            score_old = nw.NW(hmat, lecs, ':'.join(i), ':'.join(j))[2]
            score_new = nwx.NW(hmat, decs, i, j)[2]
            assert abs(score_old - score_new) <= 1e-2
            score_old_loc = nw.NW(hmat, lecs, ':'.join(i), ':'.join(j),
                                  local=True, localize=True)[2]
            score_new_loc = nwx.NW(hmat, decs, i, j, localize=True)[2]
            assert abs(score_old_loc - score_new_loc) <= 1e-2


####################
# Testing wrappers #
####################

def test_alignESS_wrap():
    """Test the correct function of the
     alignment wrapper."""
    s1 = '1.2.1:2.7.2:5.4.2:4.2.1:2.7.1'
    s2 = '9.9.9:3.2.1:5.3.1:3.1.3:4.1.2:1.2.1:5.4.2:4.2.1:2.7.1' + \
         ':1.2.7:2.3.1:1.8.1'
    ss1, ss2, s = alignESS.alignESS(s1, s2)
    ls1, ls2, ls = alignESS.alignESS(s1, s2, localize=True)
    assert ls < s
    assert len(ls2) < len(ss1)


#################################
# Testing multiple comparissons #
#################################

@skipifnotfiles
def test_indvsalldb_len(seqs, decs, hmat):
    """Test for the correct length of ind_vs_alldb

    """
    subs = seqs[:100]
    scores0 = nwx.ind_vs_alldb(0, subs, hmat, decs)
    aligns0 = nwx.ind_vs_alldb(0, subs, hmat, decs, oscore=False)
    scores89 = nwx.ind_vs_alldb(89, subs, hmat, decs)
    scoreswhole = nwx.ind_vs_alldb(89, subs, hmat, decs, wholedb=True)
    assert len(scores0) == 1
    assert len(scores89) == 1
    assert len(scoreswhole) == 1
    assert len(scores0[0]) == 99
    assert len(scores89[89]) == 10
    assert len(scoreswhole[89]) == 100
    # oscore = False result length
    assert len(aligns0) == len(scores0)
    aln = list(list(aligns0.values())[0].values())[0]
    assert len(aln) == 3
    assert type(aln[0]) is str
    assert type(aln[2]) is float


@skipifnotfiles
def test_seqvsdb_len(seqs, decs, hmat):
    """Test for the correct length output in seqs_vs_db output"""
    subs = seqs[:100]
    n = len(subs)
    seq1 = seqs[50]
    scores = nwx.seq_vs_db(seq1, subs, hmat, decs)
    aligns = nwx.seq_vs_db(seq1, subs, hmat, decs, oscore=False)
    assert len(scores) == n
    sco_thres = nwx.seq_vs_db(seq1, subs, hmat, decs, thres=0.3)
    assert len(sco_thres) < n
    # oscore = False result lenght
    assert len(scores) == len(aligns)
    aln = list(aligns.values())[0]
    assert len(aln) == 3
    assert type(aln[0]) is str
    assert type(aln[2]) is float


@skipifnotfiles
def test_dbvsdb_len(seqs, decs, hmat):
    """Test for the correct length in db_vs_db output
    """
    subs1 = seqs[:50]
    n1 = len(subs1)
    subs2 = seqs[20:60]
    n2 = len(subs2)
    scores = nwx.db_vs_db(subs1, subs2, hmat, decs)
    aligns = nwx.db_vs_db(subs1, subs2, hmat, decs,
                          oscore=False)
    assert len(scores) == n1
    assert len(scores[0]) == n2
    # oscore=False
    assert len(scores) == len(aligns)
    assert len(scores[0]) == len(aligns[0])
    # threshold
    sco_thres = nwx.db_vs_db(subs1, subs2, hmat, decs,
                             thres=0.3)
    assert len(sco_thres) < n1
    i = list(sco_thres.keys())[0]
    assert len(sco_thres[i]) < n2


@skipifnotfiles
def test_listvsalldb_len(seqs, decs, hmat):
    """Test for the correct length in list_vs_alldb output
    """
    subs = seqs[:100]
    indices = [10, 5, 80]
    scores = nwx.list_vs_alldb(indices, subs, hmat, decs)
    assert len(scores) == 3
    assert len(scores[10]) == 89
    assert len(scores[80]) == 19

    scoreswhole = nwx.list_vs_alldb(indices, subs, hmat, decs, wholedb=True)
    assert len(scoreswhole) == 3
    assert len(scoreswhole[5]) == 100
    assert len(scoreswhole[80]) == 100

    scores = nwx.list_vs_alldb(indices, subs, hmat, decs, thres=0.3)
    assert len(scores) == 3
    assert len(scores[10]) < 89
    assert len(scores[80]) < 19
    aligns = nwx.list_vs_alldb(indices, subs, hmat, decs, thres=0.3,
                               oscore=False)
    assert len(aligns) == 3
    assert len(aligns[10]) < 89
    assert len(aligns[80]) < 19

    scoreswhole = nwx.list_vs_alldb(indices, subs, hmat, decs, wholedb=True,
                                    thres=0.3)
    assert len(scoreswhole) == 3
    assert len(scoreswhole[5]) < 100
    assert len(scoreswhole[80]) < 100


@skipifnotfiles
def test_localize_len(seqs, decs, hmat):
    """Test that the option localize generates more results
    than the normal option"""
    # alldb comp
    r1 = nwx.alldb_comp(seqs, hmat, decs, thres=0.3)
    r2 = nwx.alldb_comp(seqs, hmat, decs, thres=0.3,
                        localize=True)
    s1 = [s for sub in r1.values() for s in sub.values()]
    s2 = [s for sub in r2.values() for s in sub.values()]
    n1 = len(s1)
    n2 = len(s2)
    assert n1 < n2
    assert sum(s1) / n1 > sum(s2) / n2

###############
# Testing I/O #
###############


@skipifnotfiles
def test_loadtxt():
    fname1 = 'tests/random_frag_ind.txt'
    fname2 = 'tests/random_frag.txt'

    ess = alignESS.load_txt(fname1, index=True)
    assert type(ess) == tuple
    assert len(ess) == 2
    es = ess[0][0]
    assert type(es) == list
    ess = alignESS.load_txt(fname2, index=False)
    assert type(ess) == list


@skipifnotfiles
def test_loadsqlite():
    sqlname = 'tests/nr_part.db'
    ess = alignESS.load_sqlite(sqlname)
    assert type(ess) == tuple
    assert len(ess) == 2
    assert type(ess[0][0]) == list


def test_storedict_empty():
    "test for store_dict empty"
    scores = {}
    fname = 'pru.txt'
    nwx.store_dict(fname, scores)


#####################
# Testing utilities #
#####################


@skipifnotfiles
def test_esstype():
    ess_type = alignESS._ess_type
    assert ess_type('tests/random_frag.txt') == 'text_noind'
    assert ess_type('tests/random_frag_ind.txt') == 'text_ind'
    assert ess_type('tests/part_name.txt') == 'text_ind'
    assert ess_type('tests/nr_part.db') == 'sqlite'
    assert ess_type('tests/test_nwx.py') == 'file'
    assert ess_type('1.2.3:1.3.4:3.5.-') == 'ess'


@xfail(reason='Must raise exception', strict=True)
def test_esstype_fail():
    ess_type = alignESS._ess_type
    assert ess_type('1.2.3')


@skipif(not os.path.exists('tests/multi.txt'),
        reason='File not found')
def test_loadmulti():
    load_multi = alignESS._load_multi
    ess, names = load_multi('tests/multi.txt')
    assert type(ess[0]) == list
    assert type(names[0]) == str
    assert len(ess) == 8
    assert len(names) == 8
    assert len(ess[0]) == 6
    assert 'Glycolysis' in names[3]


################################
# Testing command line options #
################################

def test_localize_parser():
    """Test that the localize option returns True or False correctly"""
    s1 = '1.2.1:2.7.2:5.4.2:4.2.1:2.7.1'
    s2 = '9.9.9:3.2.1:5.3.1:3.1.3:4.1.2:1.2.1:5.4.2:4.2.1:2.7.1' + \
         ':1.2.7:2.3.1:1.8.1'
    parser = alignESS.arg_parser(get_parser=True)
    # pair command
    args = parser.parse_args(['pair', s1, s2])
    assert args.localize is False
    args = parser.parse_args(['pair', s1, s2, '-l'])
    assert args.localize is True
    alignESS.main_pair(args)
    # dbalign
    args = parser.parse_args(['dbalign', 'tests/nr_part.db'])
    assert args.localize is False
    args = parser.parse_args(['dbalign', 'tests/nr_part.db', '-l'])
    assert args.localize is True
    alignESS.main_db(args)


def test_dbalign_align_option():
    """Test that the align option returns True or False correctly"""
    parser = alignESS.arg_parser(get_parser=True)
    # dbalign
    args = parser.parse_args(['dbalign', 'tests/nr_part.db'])
    assert args.align is False
    args = parser.parse_args(['dbalign', 'tests/nr_part.db', '-align'])
    assert args.align is True
    alignESS.main_db(args)


@skipifnotfiles
def test_singledb():
    """Using only one database yielded an error because the program tries
    to use the value in args.db2 variable. In this case this value is None"""
    parser = alignESS.arg_parser(get_parser=True)
    args = parser.parse_args(['dbalign', 'tests/nr_part.db'])
    alignESS.main_db(args)
    args = parser.parse_args(['dbalign', 'tests/part_name.txt'])
    alignESS.main_db(args)
    # comment if the test become slow
    args = parser.parse_args(['dbalign', 'tests/random_frag.txt'])
    alignESS.main_db(args)
    args = parser.parse_args(['dbalign', 'tests/random_frag_ind.txt'])
    alignESS.main_db(args)


@skipifnotfiles
def test_doubledb():
    """Using 2 databases should not give any error"""
    parser = alignESS.arg_parser(get_parser=True)
    args = parser.parse_args(['dbalign', 'tests/nr_part.db', '-db2',
                              'tests/nr_part.db'])
    alignESS.main_db(args)
    args = parser.parse_args(['dbalign', 'tests/part_name.txt', '-db2',
                              'tests/nr_part.db'])
    alignESS.main_db(args)
    args = parser.parse_args(['dbalign', 'tests/part_name.txt', '-db2',
                              'tests/random_frag.txt'])
    alignESS.main_db(args)
