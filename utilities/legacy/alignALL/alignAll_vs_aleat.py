#!/usr/bin/python
# 
#
# @uthor: acph
#  
"""Script to compare all the database vs the database of ESS versus an shuffle generated ESS database

"""

from multiprocessing import Pool
import sqlite3 as s3
from numpy import load
from nw_ec_aling import *



def load_aleat_seqs(aleat_seqs_file):
    """Loads a list of enzimatic step sequences.
    
    Arguments:
    - `aleat_seqs_file`: A text file containning a list of enzymatic step sequences in simple format. Each line contains one sequence.
    """
    aleat_seqs = []
    with open(aleat_seqs_file, 'r') as inf:
        for line in inf:
            seq = line.strip()
            aleat_seqs.append(seq)

    return aleat_seqs



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



def one2all_aleat(seq):
    """Creates the pair alingments of seq, versus the global list aleat seqs
    
    Arguments:
    - `seq`: an enzymatic step sequence, seq is a tuple (seqid, seq)
    """
    id1 = seq[0]
    seq1 = seq[1]

    results = [0] * len(aleat_seqs)

    outf = open('{}.txt'.format(id1), 'w')
    
    print 'Prossesing sequence id = {}'.format(id1)
    for i in xrange(len(aleat_seqs)):
        seq2 = aleat_seqs[i]
        s1, s2, score = NW(matrix, ecs, seq1, seq2)
        line = '{:.5}\n'.format(score)
        results[i] = line

    results = ''.join(results)
    outf.write(results)
    outf.close()


                    

if __name__ == '__main__':
    
    """ Usage

$ python alignAll_vs_aleat.py seq_database aleat_seqs_file substitution_matrix proceses_pool

    """

    from sys import argv
    from time import time
    from numpy import load
    
    t1 = time()
    
    db = s3.connect(argv[1]) # Database
#    seqs = retrieve_seq(db, 'seqs',  where="where sp='eco'")
    seqs = retrieve_seq(db, 'nrseqs')
    aleat_seqs = load_aleat_seqs(argv[2])
    
    db.close()
    # substitution matrix
    m = load(argv[3])
    matrix = m['matrix']
    ecs = m['ecs']
    ecs = list(ecs)
    # processes
    p = int(argv[4])

    pool = Pool(processes=p)
    pool.map(one2all_aleat, seqs)

    

#    map(one2all, seqs)

    elapsed = time() - t1
    elapsed = elapsed / 60
    print "Elapsed: {} minutes".format(elapsed)





