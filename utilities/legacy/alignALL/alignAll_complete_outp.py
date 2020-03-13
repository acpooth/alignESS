#!/usr/bin/python
# 
#
# @uthor: acph
#  
"""Script to do all vs all pairwise alignment of a list of EC sequences obtainde from a sqlite db. Parallel version
"""

from multiprocessing import Pool
import sqlite3 as s3
from numpy import load
from nw_ec_aling import *


def retrieve_seq(db,table, seqs='ec3', where=''):
    """Retrive all sequences from the sqli
te3 seqs database file. Returns a tuple of tuples (seq-id, seq)

db = sqlite3 database object 
table = Specifies the table to retrieve the sequences
cepted values are: seqs = Specifies the sequences type to retrieve. ec3(default), ec4, gene, arch
where = adds a sqlite where statement to retrive the sequences

retrieve_seq(args) => ((seq_id, seq),) 
"""
    if table == 'seqs':
        q = "id, " + seqs
    elif table == 'nrseqs':
        q = "nrid, " + seqs

    query = db.execute("SELECT {q} FROM {table} {where}".format(q=q, table=table, where=where))
    return tuple(query.fetchall())



def one2all(seq):
    """Creates the pair alingments of seq, versus the global list seqs
    
    Arguments:
    - `seq`: an enzymatic step sequence, seq is a tuple (seqid, seq)
    """
    index = seqs.index(seq) + 1 # index in seqs. seqs = list of al
# sequences to run
    rang = range(index)
    results = [0] * index
    id1 = seq[0]
    seq1 = seq[1]
    
    outf = open('{}.txt'.format(id1), 'w')

    print 'Prossesing sequence id = {}'.format(id1)
    for i  in rang:
        seq2 = seqs[i][1]
        id2 = seqs[i][0]
        s1, s2, score = NW(matrix, ecs, seq1, seq2)
        line = '{}\t{}\t{:.5}\t{}\t{}\n'.format(id1, id2, score, s1, s2)
        results[i] = line

    results = ''.join(results)
#    return results
    outf.write(results)
    outf.close()




if __name__ == '__main__':
    """ Usage

$ python alignAll.py seq_database substitution_matrix proceses_pool where_statment

    """

    from sys import argv
    from time import time
    from numpy import load
    
    t1 = time()
    
    db = s3.connect(argv[1]) # Database
#    seqs = retrieve_seq(db, 'seqs',  where="where sp='eco'")
    seqs = retrieve_seq(db, 'nrseqs')
    db.close()
    print len(seqs)
    # substitution matrix
    m = load(argv[2])
    matrix = m['matrix']
    ecs = m['ecs']
    ecs = list(ecs)
    # processes
    p = int(argv[3])

    pool = Pool(processes=p)
    pool.map(one2all, seqs)


    

#    map(one2all, seqs)

    elapsed = time() - t1
    elapsed = elapsed / 60
    print "Elapsed: {} minutes".format(elapsed)

    












