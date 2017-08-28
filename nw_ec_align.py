#!/usr/bin/python
#
#
# @uthor: acph (dragopoot@gmail.com)
#

"""Program to align a pair of enzimatic step sequences  using an implementation of the Needleman-Wuncsh  dynamic programing algorithm """

import argparse
from numpy import zeros, arange, random, sign, mean, ones, load
import re


# Load matrices
fol = __file__.split('/')
if len(fol) == 1:
    fol = '.'
else:
    fol = '/'.join(fol[:-1])
#mat = load(fol + '/entropy_matrix.npz')
#mat = mat['matrix']

hmat = load(fol + '/entropy_matrix_hierarchy.npz')
ecs = hmat['ecs']
hmat = hmat['matrix']


# -- End load matrices

def scoring(mat, alp, seq1, seq2, gap=1, fhomo=0.95, fpengap=0.05, local=False):
    """ This function evaluates an enzymatic step sequences pair alignment. The sequences need to be aligned. This function tries to recall the fitness function in Ortegon-Cano, 2011. 

    Arguments:
    - `mat`: score matrix
    - `alp`: alphabet, index of score matrix
    - `s1`: enzimatic step sequence 1, ':' delimited
    - `s2`: enzimatic step sequence 2, ':' delimited
    - `gap`: gap penalization
    - `fhomo`: Weight factor to multiplicate the homogeneity fraction of score
    - `fpengap`: Weight factor to multiplicate the gap penalization
    - `local`: If True, evaluates the score in a local fassion, rather than global, i.e.  eliminates from the alignment the initial and final gaps of the shortest sequence and the respective positions in the largest sequence.
    """
    s1 = seq1.strip(' :\n\r')
    s2 = seq2.strip(' :\n\r')
    s1 = s1.split(":")
    s2 = s2.split(":")
    assert len(s1) == len(
        s2), 'The sequences aligned must be of the same length'

    # -- Mean column homogeneity evaluation
    if local:
        length = len(s1)
        ii = 0
        for i in range(length):
            if s1[i] != '-.-.-' and s2[i] != '-.-.-':
                ii = i    # initial index
                break
        fi = length
        for i in range(length)[::-1]:
            if s1[i] != '-.-.-' and s2[i] != '-.-.-':
                break
            fi = i      # final index
        # trimming the sequences
        s1 = s1[ii:fi]
        s2 = s2[ii:fi]

    length = len(s1)  # Alignment length
    homo_col = []
    for i in range(length):
        if s1[i] == '-.-.-' or s2[i] == '-.-.-':  # and s1[i] != s2[1]
            homo_col.append(gap)
        elif s1[i] == '...' or s2[i] == '...':
            pass
        elif s1[i] == '' or s2[i] == '':
            pass
        else:
            x = alp.index(s1[i])
            y = alp.index(s2[i])
            homo_col.append(mat[x, y])  # pondered by 0.6 as in

    homo = mean(homo_col)
    # -- Column increment penalization > here isnt necesary,isn't it?
    # ls1 = len([s for s in s1 if s != '-.-.-'])
    # ls2 = len([s for s in s2 if s != '-.-.-'])
    # incpen = (length - max(ls1, ls2))/10.
    # -- Gap penalization
    seq1 = seq1.strip('.-:')  # remove initial and final gaps
    seq2 = seq2.strip('.-:')
    gap = re.compile(r"(:-.-.-)+")  # gap regex

    gaps = 0  # gap penalty value
    for s in (seq1, seq2):
        sgaps = []  # list of gap blocks
        for g in gap.finditer(s):
            g = g.group()
            g_lenght = g.count("-.-.-")
            sgaps.append(g_lenght)
        if sgaps == []:
            gaps += 0
        else:
            gaps += float(len(sgaps)) / sum(sgaps)  # num of gap blocks /
# num of independent blocks per sequence
    pengap = gaps / 2  # mean , between number of sequences

    # Final score calculation
    score = (homo * fhomo) + (pengap * fpengap)  # + incpen
    return score


def diagonal_list(l1, l2):
    """Creates a list of diagonals in a matrix. Index begin = 1. Used to facilitate the calculation of the scoring and arrow matrices.

    Arguments:
    - `l1`: lenght of sequence 1
    - `l2`: lenght of sequence 2
    """
    dlist = []
    for i in range(l1 + l2 - 1):
        st1 = min(i + 1, l1)
        sp1 = max(1, i - l2 + 2)
        st2 = max(1, i - l1 + 2)
        sp2 = min(i + 1, l2)
#        print st1,sp1,st2,sp2
        dlist.append((arange(st1, sp1 - 1, -1), arange(st2, sp2 + 1)))
    return dlist


def FastSubValues(mat, alp, s1, s2):
    """Creates a matrix with the same dimentions of scoring matrix and arrow matrix that contains the precomputed values of the alignment of each comparation in both sequences.

    Arguments:
    - `mat`: substitution matrix
    - `alp`: alphabeth list, index of the subvalue matrix
    - `s1`: enzimatic step sequence 1, ':' delimited
    - `s2`: enzimatic step sequence 2, ':' delimited
    """
    s1 = s1.strip(': \n\r')
    s2 = s2.strip(': \n\r')
    s1 = s1.split(':')
    s2 = s2.split(':')
#    alp = tuple(alp)
    l1, l2 = len(s1), len(s2)
    subvals = zeros((l1 + 1, l2 + 1))  # substitution values matrix
    # Convert the sequences to sequences of index
    si1 = zeros(l1, int)
    si2 = zeros(l2, int)
    for i in range(l1):
        si1[i] = alp.index(s1[i])
    for i in range(l2):
        si2[i] = alp.index(s2[i])
    # Fills subvalue matrix, this initiates in the position (1,1)
    # similar, to scoring and arrow matrices
    for i in range(1, l1 + 1):
        subvals[i, 1:] = mat[[si1[i - 1]] * l2, si2]
    return subvals


def FastNW(subvals, s1, s2, gap=0.9):
    """Perform dynamic programing alignment of EC numbers sequences. Creates the scoring and arrow matrices using subvals matrix and diagonal arrays.

    Returns the dynamic programing matrix and the arrow matrix used to trace back the alignment.

    Arguments:
    - `subvals`: matrix of precomputed values for the alignment of each pair of characters
    - `s1`: enzymatic step sequence 1, ':' delimited
    - `s2`: enzymatic step sequence 2, ':' delimited
    - `gap`: gap penalties, default = 1
    """
    s1 = s1.strip(': \n\r')
    s2 = s2.strip(': \n\r')
    s1 = s1.split(':')
    s2 = s2.split(':')
    l1, l2 = len(s1), len(s2)
    # Create the score and aroow matrices
    scoremat = zeros((l1 + 1, l2 + 1))
    arrow = zeros((l1 + 1, l2 + 1))
    # Create first row and first column with gaps
    scoremat[0] = arange(l2 + 1) * gap
    scoremat[:, 0] = arange(l1 + 1) * gap
    arrow[0] = ones(l2 + 1)
    # Compute diagonal list
    dlist = diagonal_list(l1, l2)
    # fill the matrix
    for i in dlist:
        li = len(i[0])
        f = zeros((3, li))  # results of the tree posibles xhoices
# for each  value in the diagonal
        x, y = i[0], i[1]
        f[0] = scoremat[x - 1, y] + gap
        f[1] = scoremat[x, y - 1] + gap
        f[2] = scoremat[x - 1, y - 1] + subvals[i]
        f -= 0.001 * sign(f) * random.ranf(f.shape)  # for randomly
# select from a tie
        mini = f.min(0)
        minpos = f.argmin(0)
        scoremat[i] = mini
        arrow[i] = minpos
    return scoremat, arrow, mini


def backtrace(arrow, s1, s2):
    """Reads the arrow matrix and return the aligned sequences. Global, N-W

    Arguments:
    - `arrow`: arrow matrix generated by dynamic programing
    - `s1`: enzymatic step sequence, ':' delimited
    - `s2`: enzymatic step sequence, ':' delimited
    """
    # Transform the sequences strings into EC numbers lists
    s1 = s1.strip(': \n\r')
    s2 = s2.strip(': \n\r')
    s1 = s1.split(':')
    s2 = s2.split(':')
    st1, st2 = [], []  # aligned sequences
    ok = True
    v, h = arrow.shape
    v -= 1
    h -= 1
    while ok:
        if arrow[v, h] == 0:  # vertical best result, s1
            st1.append(s1[v - 1])
            st2.append('-.-.-')
            v -= 1
        elif arrow[v, h] == 1:  # horizontal best result, s2
            st1.append('-.-.-')
            st2.append(s2[h - 1])
            h -= 1
        elif arrow[v, h] == 2:  # diagonal best result, s1,s2 aligned
            st1.append(s1[v - 1])
            st2.append(s2[h - 1])
            v -= 1
            h -= 1
        if v == 0 and h == 0:
            ok = False
    # reverse sequences
    st1 = ':'.join(st1[::-1])
    st2 = ':'.join(st2[::-1])
    return st1, st2


def NW(mat, ecs, s1, s2, gap=0.9, local=False, localize=False, nws=False):
    """
    NW alignment function. Creates a pairwise alignment using a Needelman-Wunsh algorithm.

    Arguments:
    - `mat`: EC number (3 levels) substitution matrix
    - `ecs`: List of ec numbers that represent the labels of the matrix
    - `s1`: EC numbers sequence 1, : delimited
    - `s2`: EC numbers sequence 2, : delimited
    - `gap`: gap penalty
    - `local`: if local = True, then the score of the alignment is calculated 
               localy, i.e. only en the part of the alignment covered by the 
               shortest sequence
    - `localize`: if localize = True, then, the function returns only the 
               fragment of the alignment covered by the shortest sequence
    - `nws`: if true, the function returns the scores as is returned by the 
               NW algorithm
    """
    submat = FastSubValues(mat, ecs, s1, s2)  # create sub score matrix
    sco, arr, mini = FastNW(submat, s1, s2, gap=gap)  # create score and
# arrow matrices
    s1, s2 = backtrace(arr, s1, s2)  # backtrace alingment
    # scoring : the scoring is doing using the entropy schema to
    # normalize the score in a range 0 to 1 and to generate similar
    # results to the genetic algorithm
    if nws:                    # NW score
        if localize:
            s1 = s1.split(':')
            s2 = s2.split(':')
            length = len(s1)
            ii = 0
            for i in range(length):
                if s1[i] != '-.-.-' and s2[i] != '-.-.-':
                    ii = i    # initial index
                    break
            fi = length
            for i in range(length)[::-1]:
                if s1[i] != '-.-.-' and s2[i] != '-.-.-':
                    break
                fi = i      # final index

            s1 = s1[ii:fi]
            s2 = s2[ii:fi]
            s1 = ':'.join(s1)
            s2 = ':'.join(s2)
        return s1, s2, mini[0]
    scoring_gap = 1
    score = scoring(mat, ecs, s1, s2, gap=scoring_gap, local=local)
    if localize:
        s1 = s1.split(':')
        s2 = s2.split(':')
        length = len(s1)
        ii = 0
        for i in range(length):
            if s1[i] != '-.-.-' and s2[i] != '-.-.-':
                ii = i    # initial index
                break
        fi = length
        for i in range(length)[::-1]:
            if s1[i] != '-.-.-' and s2[i] != '-.-.-':
                break
            fi = i      # final index

        s1 = s1[ii:fi]
        s2 = s2[ii:fi]
        s1 = ':'.join(s1)
        s2 = ':'.join(s2)
    return s1, s2, score


################
# # -          #
# # --- legacy #
# # -          #
################

def scoring_l(mat, alp, seq1, seq2, gap=1, fhomo=0.95, fpengap=0.05):
    """ This function evaluates an enzymatic step sequences pair alignment. The sequences need to be aligned. This function tries to recall the fitness function in Ortegon-Cano, 2011. LEGACY

    Arguments:
    - `mat`: score matrix
    - `alp`: alphabeth, index of score matrix
    - `s1`: enzimatic step sequence 1, ':' delimited
    - `s2`: enzimatic step sequence 2, ':' delimited
    - `gap`: gap penalization
    - `fhomo`: Weight factor to multiplicate the homogeneity fraction of score
    - `fpengap`: Weight factor to multiplicate the gap penalization
    """
    s1 = seq1.strip(' :\n')
    s2 = seq2.strip(' :\n')
    s1 = s1.split(":")
    s2 = s2.split(":")
    assert len(s1) == len(
        s2), 'The sequences aligned must be of the same length'

    length = len(s1)  # Alignment length

    # -- Mean column homogeneity evaluation
    homo_col = []
    for i in range(length):
        if s1[i] == '-.-.-' or s2[i] == '-.-.-':  # and s1[i] != s2[1]:
            homo_col.append(gap)
        elif s1[i] == '...' or s2[i] == '...':
            pass
        elif s1[i] == '' or s2[i] == '':
            pass
        else:
            x = alp.index(s1[i])
            y = alp.index(s2[i])
            homo_col.append(mat[x, y])  # pondered by 0.6 as in
# Ortegon-Cano 2011.This value is taken from the multiple alignment AG
    homo = mean(homo_col)

    # -- Column increment penalization > here isnt necesary

    # -- Gap penalization
    seq1 = seq1.strip('.-:')  # remove initial and final gaps
    seq2 = seq2.strip('.-:')
    gap = re.compile(r"(:-.-.-)+")  # gap regex

    # gaps = [] # count of gaps blocks,
    # # lenght = number of blocks in alignment
    # # sum = number of individual gaps in the alignment
    # for seq in (seq1,seq2):
    #     for g in gap.finditer(seq):
    #         g = g.group()
    #         g_lenght = g.count("-.-.-")
    #         gaps.append(g_lenght)
    # total_gaps = sum(gaps)
    # total_blocks = len(gaps)
    # pengap = float(total_gaps)/total_blocks/total_gaps
    # #to minimize
    # pengap = 1 - pengap

    # alt
    gaps = 0  # gap penalty value
    for s in (seq1, seq2):
        sgaps = []  # list of gap blocks
        for g in gap.finditer(s):
            g = g.group()
            g_lenght = g.count("-.-.-")
            sgaps.append(g_lenght)
        if sgaps == []:
            gaps += 0
        else:
            gaps += float(len(sgaps)) / sum(sgaps)  # num of gap blocks /
# num of independent blocks per sequence
    pengap = gaps / 2  # mean

    # Final score calculation

    score = (homo * fhomo) + (pengap * fpengap)
    return score


def NW_l(mat, ecs, s1, s2, gap=0.9):
    """
    NW alignment function. Creates a pairwise alignment using a Needelman-Wunsh algorithm. LEGACY

    Arguments:
    - `mat`: EC number (3 levels) substitution matrix
    - `ecs`: List of ec numbers that represent the labels of the matrix
    - `s1`: EC numbers sequence 1, : delimited
    - `s2`: EC numbers sequence 2, : delimited
    - `gap`: gap penalty 
    """
    submat = FastSubValues(mat, ecs, s1, s2)  # create sub score matrix
    sco, arr, mini = FastNW(submat, s1, s2, gap=gap)  # create score and
# arrow matrices
    s1, s2 = backtrace(arr, s1, s2)  # backtrace alingment
    # scoring : the scoring is doing using the entropy schema to
    # normalize the score in a range 0 to 1 and to generate similar
    # results to the genetic algorithm
    scoring_gap = 1
    score = scoring(mat, ecs, s1, s2, gap=scoring_gap)
    return s1, s2, score  # , mini # uncoment mini to show the score
# frome the scoring matrix


def arg_parser():
    parser = argparse.ArgumentParser(
        description='Needleman-Wuncsh algorithm for alignment of ESS (Enzymatic Step Sequence)')
    parser.add_argument(
        'ess1', type=str, help='Sequence of 3 levels EC numbers. Colon separated. (1.2.3:3.5.-:...:9.9.9)')
    parser.add_argument(
        'ess2', type=str, help='Sequence of 3 levels EC numbers. Colon separated. (1.2.3:3.5.-:...:9.9.9)')
    parser.add_argument('--gap', type=float, default=0.9,
                        help='Gap penalization (from 0 to 1). Default = 0.9')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    from sys import stdout

    args = arg_parser()
    ecs = list(ecs)
    result = NW(hmat, ecs, args.ess1, args.ess2, gap=args.gap)

    stdout.write('ess1:\t{}\n'.format(result[0]))
    stdout.write('ess2:\t{}\n'.format(result[1]))
    stdout.write('score = {}\n'.format(result[2]))
