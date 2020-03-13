#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       @uthor: acph
# 02.09.2011 08:48:33

import subprocess
from time import time	
import sqlite3 as s3
from multiprocessing import Pool



def pair_alignment(seq1, seq2):
    """Generates the aligment using GA of the seq1 and seq2 enzymatic steps sequences

    pair_alignment(seq1, seq2) => list[seq1-al, seq2-al, fitness] 
"""
    cline = "./AlineaPares 100 100 0.9 0.01 0.05 0.95 1 -5 .25".split()
    cline += [seq1, seq2]
    child = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    handler = child.stdout
    lines = handler.read().split('\n')
    return lines[-4:-1]


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


 
def seq_vs_all(seq):
    """Do pairwise alignments of the sequence tuple (seq =(id, seq)) versus all the sequences in the tuple of seqeunces (tuple_seqs) generated by retrive_seqs*=()

Returns => a set of strings. Each string contains the infromation seq1, seq2 and fitness tab delimited.
"""
    print 'Prossesing sequence id = {} ...'.format(seq[0])
    ind = seqs.index(seq)
    rang = range(ind)
    results = [''] * (ind)
    for i in rang:
        fit = pair_alignment(seq[1], seqs[i][1])
        fit = fit[-1].split()[1]
        results[i] = "{0}\t{1}\t{2}".format(seq[0], seqs[i][0], fit)
    with open('{}.txt'.format(seq[0]), 'w') as outf :
        outf.write('\n'.join(results))
    return 1
    



if __name__ == '__main__':
    '''main for compare all sequences from de database
    argument1 = seqs db file (seqs.db)
    argument2 = seqid init
    argument3 = seqid final
    
needs to be in the same directory as AlineaPares exe '''
    from sys import argv

# data
    t1 = time()
    db = s3.connect(argv[1])
    try:
        where = 'WHERE nrid>{}'.format(argv[2])
    except:
        print "No inital value especified ..."
        print "Use '0' to start from the begining"
        print "Exit!"
        exit()

    try:
        where += ' and nrid<{}'.format(argv[3])
    except:
        print 'No end value used...  i will run until de last element in nrseqs'
    
    
    seq = retrieve_seq(db, 'nrseqs', where=where)
    seqs = retrieve_seq(db, 'nrseqs')
    db.close()


# process    
    pool = Pool(processes= 4)
    res = pool.map(seq_vs_all, seq)

# result storage
    
# print report
    elapsed = time() - t1

    print "Time elpased: {0} minutes".format(elapsed/60)

    