#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
# Name:     
# Purpose:  Script to run al the comparisson of an ESS database
# 
# @uthor:   acph - dragopoot@gmail.com
#
# Created:     
# Copyright:   (c) acph 2014
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#------------------------------
"""Script to run al the comparisson of an ESS database.
Align the database vs itself and compare with random
generated databases"""

import sys
import os
import ctypes
import sqlite3 as s3
import numpy as np
from multiprocessing import Pool, Array
from functools import partial
# pysrc = '/home/acph/pysrc'
# ESS = '/home/acph/Proyectos/ESS_AlingmentsDP'
# kegg2seq = '/home/acph/Dropbox/scripts/kegg2seq_py'
DP = '/home/acph/Dropbox/scripts/Alignment_EC_DP'
# lista = [pysrc, ESS, kegg2seq, DP]
sys.path = sys.path + [DP]
import nw_ec_align as NW
# import anubis as A
# import scores_tools as S
# import cluster_tools as C
# from kegg_maps import keggmaps

# Load files
db = s3.connect("/home/acph/KEGG2seq/kegg2seq_py-bfs/kegg2sec/Db/seqs.db")

hmat = NW.hmat
ecs = list(NW.ecs)

########################
# Load data from files #
########################
def load_random_db(filename):
    """Loads the random Db text file. Returns a dictionary: keys = id
    values = ec3 ESS
    
    Arguments:
    - `filename`: database filename
    """
    essdic = {}
    with open(filename) as inf:
        for line in inf:
            line = line.split('\t')
            _id = int(line[0])
            ess = line[1]
            essdic[_id] = ess
    return essdic

def get_real_ess(db=db):
    """Returns a dictionary. keys = nrid. values = another dictionary. keys
    ec3, maps; values ec3 ESS and a set of metabolic maps, 5 digits.
    
    Arguments:
    - `db`: sqlite3 ESS database object
    """
    handler = db.execute("SELECT nrid, ec3, maps, len FROM nrseqs")
    nrdic = {}
    for query in handler:
        nrid = query[0]
        ec3 = query[1]
        maps = query[2].split()
        maps = set([m[3:] for m in maps])
        _len = query[3]
        nrdic[nrid] ={'ec3': ec3,
                      'maps': maps,
                      'len': _len}
    return nrdic

# - Constants database file
realdic = get_real_ess()

#########################
# -                     #
# -- Align real Vs real #
# -                     #
#########################

def align(ess_pair, nrdic):
    """Wraper for easy align ESS
    
    Arguments:
    - `ess_pair`: tuple. Pair of ess ids
    """
    assert len(ess_pair) == 2, "Must be a pair of ess id"
    id1, id2 = ess_pair
    ess1 = nrdic[id1]
    ess2 = nrdic[id2]
    result = NW.NW(hmat, ecs, ess1, ess2, local=True)
    score = result[2]
    data = np.array([id1, id2, score], dtype = np.float32)
    return data

def fill_mat(ess_pair,  nrdic):
    """
    """
    id1, id2, score = align(ess_pair, nrdic)
    scomat[id1-1, id2-1] = score
    
def align_real():
    """Aling all the sequenes in the sqlite database
    
    """
    # simplify the realdic dictionary into nrdic k=nrid, v=ec3 ESS
    nrdic = {k: v['ec3'] for k, v in realdic.iteritems()}
    total = np.max(nrdic.keys())
    i1, i2 = np.triu_indices(total)
    i1 = np.int16(i1 + 1)
    i2 = np.int16(i2 + 1)
    indices = np.vstack((i1, i2)).T
    # scomat = np.zeros((total, total), dtype = np.float16)
    # alternative scomat
    shared_sco_base = Array(ctypes.c_float, total*total)
    # scomat = np.ctypeslib.as_array(shared_sco_base.get_obj())
    scomat_ = np.frombuffer(shared_sco_base.get_obj(), dtype=np.float32)
    scomat = scomat_.reshape((total, total))
    global scomat
    align_func = partial(fill_mat,  nrdic=nrdic)
    pool = Pool(processes=6)
    pool.map(align_func, indices, chunksize=1000)
    # ids1 and ids2 are not strictly necessary because data is stored in
    # indices
    ids1 = np.zeros(len(indices))
    ids2 = np.zeros(len(indices))
    scores = np.zeros(len(indices))
    for i in xrange(len(indices)):
        id1, id2 = indices[i]
        sco = scomat[id1-1, id2-1]
        ids1[i] = id1
        ids2[i] = id2
        scores[i] = sco
    out_filename = 'scores_local_scoring'
    np.savez(out_filename, ids1=ids1, ids2=ids2, scores=scores)
    # return scomat

###########################
# -                       #
# -- Align real Vs Random #
# -                       #
###########################

def align_ess_vs_real(rid, nrdic):
    """Wraper for easy align ESS
    One random ESS versus all the real
    
    Arguments:
    - `rid`: nrdic, ress id (randome ess id)
    """
    assert rid in nrdic, "ress must be in nrdic"
    ress = nrdic[rid]           # ress = random ess
    scores = np.zeros(len(realdic), dtype=np.float32)
    for nrid, v in realdic.iteritems():
        ess = v['ec3']
        result = NW.NW(hmat, ecs, ess, ress, local=True)
        score = result[2]
        scores[nrid-1] = score
    # fill mat
    scomat[:, rid-1] = scores
    
def align_allreal_vs_allrandom(rd_filename, out_filename):
    """Aling all the sequences in the database
    
    Arguments:
    - `rd_filename`: random database filename
    - `out_filename`: output file name -
                   output is a .npz with 3 arrays
                   ids1 (int16), ids2 (int16) , scores(float32)
    """
    nrdic = load_random_db(rd_filename)
    id_list = sorted(nrdic.keys())
    assert len(realdic) == len(nrdic), "databases must be same length"
    total = np.max(nrdic.keys())
    shared_sco_base = Array(ctypes.c_float, total*total)
    scomat_ = np.frombuffer(shared_sco_base.get_obj(), dtype=np.float32)
    scomat = scomat_.reshape((total, total))
    global scomat
    align_func = partial(align_ess_vs_real,  nrdic=nrdic)
    pool = Pool(processes=6)
    pool.map(align_func, id_list)

    np.save(out_filename, scomat)
    # return scomat

def main_random_vs_real():
    """Aligns the real database vs the 10 random databases
    """
    ifolder = 'databases/'
    ofolder = 'real_vs_random/'
    if not os.path.exists(ofolder):
        os.mkdir(ofolder)
    if_list = os.listdir(ifolder)
    if_list = [_if for _if in if_list if '.txt' in _if]
    if_list.sort()
    for infile in if_list:
        print "working {}...".format(infile)
        ifname = ifolder + infile
        number = infile.split('_')[2]
        oname = '{}_realVSran_scomat'.format(number)
        ofname = ofolder + oname
        align_allreal_vs_allrandom(ifname, ofname)

#######################
# -                   #
# -- Random Vs random #
# -                   #
#######################
    
def align_allran_vs_allran(rd_filename, out_filename):
    """Aling all the sequences in the database
    
    Arguments:
    - `rd_filename`: random database filename
    - `out_filename`: output file name -
                   output is a .npz with 3 arrays
                   ids1 (int16), ids2 (int16) , scores(float32)
    """
    nrdic = load_random_db(rd_filename)
    total = np.max(nrdic.keys())
    i1, i2 = np.triu_indices(total)
    i1 = np.int16(i1 + 1)
    i2 = np.int16(i2 + 1)
    indices = np.vstack((i1, i2)).T
    # scomat = np.zeros((total, total), dtype = np.float16)
    # alternative scomat
    shared_sco_base = Array(ctypes.c_float, total*total)
    # scomat = np.ctypeslib.as_array(shared_sco_base.get_obj())
    scomat_ = np.frombuffer(shared_sco_base.get_obj(), dtype=np.float32)
    scomat = scomat_.reshape((total, total))
    global scomat
    align_func = partial(fill_mat,  nrdic=nrdic)
    pool = Pool(processes=4)
    pool.map(align_func, indices, chunksize=1000)

    ids1 = np.zeros(len(indices))
    ids2 = np.zeros(len(indices))
    scores = np.zeros(len(indices))
    for i in xrange(len(indices)):
        id1, id2 = indices[i]
        sco = scomat[id1-1, id2-1]
        ids1[i] = id1
        ids2[i] = id2
        scores[i] = sco
    np.savez(out_filename, ids1=ids1, ids2=ids2, scores=scores)
    # return scomat

def main_random_vs_random():
    """Create all the files for the random databases
    random vs random
    """
    ifolder = 'databases/'
    ofolder = 'random_vs_random'
    if not os.path.exists(ifolder):
        os.mkdir(ofolder)
    if_list = os.listdir(ifolder)
    if_list = [_if for _if in if_list if '.txt' in _if]
    if_list.sort()
    for _if in if_list:
        print "working {}...".format(_if)
        ifname = ifolder + _if
        number = _if.split('_')[2]
        oname = '{}_ranVSran'.format(number)
        ofname = ofolder + oname
        align_allran_vs_allran(ifname, ofname)

        

if __name__ == '__main__':
    # - Align real vs real
    align_real()
    # - Align random vs random
    main_random_vs_random()
    # - Align random vs real
    main_random_vs_real()

    