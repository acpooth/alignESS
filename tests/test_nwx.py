#- pytest suit
#- test for nw_ec_alignx.pyx

import pytest
from random import choice, randint
import pyximport
import alignESS
import nw_ec_align as nw
pyximport.install()
import nw_ec_alignx as nwx

decs = alignESS.decs
hmat = alignESS.hmat


@pytest.fixture(scope="module")
def seqs_s():
    import sqlite3 as s3
    with s3.connect('nr.db') as db:
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


def test_indvsalldb_len(seqs, decs, hmat):
    """Test for the correct length of ind_vs_alldb

    """
    subs = seqs[:100]
    scores0 = nwx.ind_vs_alldb(0, subs, hmat, decs)
    scores89 = nwx.ind_vs_alldb(89, subs, hmat, decs)
    scoreswhole = nwx.ind_vs_alldb(89, subs, hmat, decs, wholedb=True)
    assert len(scores0) == 1
    assert len(scores89) == 1
    assert len(scoreswhole) == 1
    assert len(scores0[0]) == 99
    assert len(scores89[89]) == 10
    assert len(scoreswhole[89]) == 100


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

    scoreswhole = nwx.list_vs_alldb(indices, subs, hmat, decs, wholedb=True,
                                    thres=0.3)
    assert len(scoreswhole) == 3
    assert len(scoreswhole[5]) < 100
    assert len(scoreswhole[80]) < 100
