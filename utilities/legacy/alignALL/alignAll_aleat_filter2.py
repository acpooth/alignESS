#!/usr/bin/python
# 
#
# @uthor: acph
#  

"""FILTERED: only uses those sequences with len > 2
"""


from multiprocessing import Pool
import sqlite3 as s3
from numpy import load
from nw_ec_aling import *



def load_aleat_seqs(aleat_seqs_file, lenght=None):
    """Loads a list of enzimatic step sequences.
    
    Arguments:
    - `aleat_seqs_file`: A text file containning a list of enzymatic step sequences in simple format. Each line contains one sequence.
    - `lenght`: if not None, it takes an integer. The sequences will be filtered using only those sequences with lenght abode lenght
    """
    aleat_seqs = []
    if lenght == None:
        with open(aleat_seqs_file, 'r') as inf:
            for line in inf:
                line = line.strip()
                seq = line.split()
                aleat_seqs.append(seq)
            
        return aleat_seqs

    elif type(lenght) == int:
        with open(aleat_seqs_file, 'r') as inf:
            for line in inf:
                line = line.strip()
                if line == '': continue
                seq = line.split()
                l = len(seq[1].split(':'))
                if l <= lenght:
                    continue
                aleat_seqs.append(seq)
            
        return aleat_seqs
        
        




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
        line = '{}\t{}\t{:.5}\n'.format(id1, id2, score)
        results[i] = line

    results = ''.join(results)
#    return results
    outf.write(results)
    outf.close()




if __name__ == '__main__':
    """ Usage

$ python alignAll.py aleat_seqs substitution_matrix proceses_pool 

    """

    from sys import argv
    from time import time
    from numpy import load
    
    t1 = time()
    
    seqs = load_aleat_seqs(argv[1], 2)
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

