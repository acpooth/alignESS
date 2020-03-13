# from libcpp cimport bool
import cython
cimport cython
from cpython cimport bool
import re
import numpy as np
cimport numpy as np
from libc.stdlib cimport rand
from libc.math cimport copysignf
cdef extern from "limits.h":
    int INT_MAX

DTYPE = np.int
ctypedef np.int_t DTYPE_t
DTYPEF = np.float
ctypedef np.float_t DTYPEF_t


gapre = re.compile(r"(:-.-.-)+")  # gap regex


def diagonal_list(int l1, int l2):
    """Creates a list of diagonals in a matrix. Index begin = 1. Used to
    facilitate the calculation of the scoring and arrow matrices.

    Generator version

    Arguments:
    - `l1`: lenght of sequence 1
    - `l2`: lenght of sequence 2
    """
    cdef int st1, sp1, st2, sp2, i
    for i in range(l1 + l2 - 1):
        st1 = min(i + 1, l1)
        sp1 = max(1, i - l2 + 2)
        st2 = max(1, i - l1 + 2)
        sp2 = min(i + 1, l2)
        # print st1,sp1,st2,sp2
        yield (np.arange(st1, sp1 - 1, -1), np.arange(st2, sp2 + 1))


# Best version but range dont allow operands
def _diagonal_list(int l1, int l2):
    """Creates a list of diagonals in a matrix. Index begin = 1. Used to
    facilitate the calculation of the scoring and arrow matrices.

    Arguments:
    - `l1`: lenght of sequence 1
    - `l2`: lenght of sequence 2
    """
    dlist = []
    cdef int st1, sp1, st2, sp2, i
    for i in range(l1 + l2 - 1):
        st1 = min(i + 1, l1)
        sp1 = max(1, i - l2 + 2)
        st2 = max(1, i - l1 + 2)
        sp2 = min(i + 1, l2)
        # print st1,sp1,st2,sp2
        dlist.append((range(st1, sp1 - 1, -1), range(st2, sp2 + 1)))
    return dlist


def _FastSubValues(mat, alp, s1, s2):
    """Creates a matrix with the same dimentions of scoring matrix and arrow
    matrix that contains the precomputed values of the alignment of each
    comparation in both sequences.

    Arguments:
    - `mat`: substitution matrix
    - `alp`: alphabeth dictionary, keys are elements, values are indices
             in substitution mat.
    - `s1`: list, enzimatic step sequence 1
    - `s2`: list, enzimatic step sequence 2
    """
    cdef int l1, l2, i, si
    l1, l2 = len(s1), len(s2)
    subvals = np.zeros((l1 + 1, l2 + 1))  # substitution values matrix
    # Convert the sequences to sequences of index
    si2 = [alp[ec] for ec in s2]
    # Fills subvalue matrix, this initiates in the position (1,1)
    # similar, to scoring and arrow matrices
    for i in range(1, l1 + 1):
        si = i - 1
        subvals[i, 1:] = mat[alp[s1[si]], si2]
    return subvals


# Best version
@cython.boundscheck(False)
@cython.wraparound(False)
#@cython.nonecheck(False)
def FastSubValues(np.ndarray[DTYPEF_t, ndim=2] mat, dict alp, s1, s2):
    """Creates a matrix with the same dimentions of scoring matrix and arrow
    matrix that contains the precomputed values of the alignment of each
    comparation in both sequences.

    Arguments:
    - `mat`: substitution matrix
    - `alp`: alphabeth dictionary, keys are elements, values are indices
             in substitution mat.
    - `s1`: list, enzimatic step sequence 1
    - `s2`: list, enzimatic step sequence 2
    """
    cdef unsigned int l1, l2, i, si, j
    l1, l2 = len(s1), len(s2)
    cdef np.ndarray[DTYPE_t, ndim= 1] si2 = np.zeros(l2, DTYPE)
    cdef np.ndarray[DTYPEF_t, ndim= 2] subvals = np.zeros((l1 + 1, l2 + 1), DTYPEF)
    # Convert the sequences to sequences of index
    for j, ec in enumerate(s2):
        si2[j] = alp[ec]
    # Fills subvalue matrix, this initiates in the position (1,1)
    # similar, to scoring and arrow matrices
    for i in range(1, l1 + 1):
        si = i - 1
        subvals[i, 1:] = mat[alp[s1[si]], si2]
    return subvals


def FastNW(subvals, s1, s2, float gap=0.9):
    """Perform dynamic programing alignment of EC numbers sequences. Creates
    the scoring and arrow matrices using subvals matrix and diagonal arrays.

    Returns the dynamic programing matrix and the arrow matrix used to trace
    back the alignment.

    Arguments:
    - `subvals`: matrix of precomputed values for the alignment of each pair
                of characters
    - `s1`: list, enzymatic step sequence 1
    - `s2`: list, enzymatic step sequence 2
    - `gap`: gap penalties, default = 1
    """
    cdef int l1, l2, li
    l1, l2 = len(s1), len(s2)
    # Create the score and aroow matrices
    scoremat = np.zeros((l1 + 1, l2 + 1))
    arrow = np.zeros((l1 + 1, l2 + 1))
    # Create first row and first column with gaps
    scoremat[0] = np.arange(l2 + 1) * gap
    scoremat[:, 0] = np.arange(l1 + 1) * gap
    arrow[0] = np.ones(l2 + 1)
    # Compute diagonal list
    dlist = diagonal_list(l1, l2)
    # fill the matrix
    for i in dlist:
        li = len(i[0])
        f = np.zeros((3, li))  # results of the tree posibles xhoices
# for each  value in the diagonal
        x = i[0]
        y = i[1]
        f[0] = scoremat[x - 1, y] + gap
        f[1] = scoremat[x, y - 1] + gap
        f[2] = scoremat[x - 1, y - 1] + subvals[i]
        f -= 0.001 * np.sign(f) * np.random.ranf(f.shape)  # for randomly
# select from a tie
        mini = f.min(0)
        minpos = f.argmin(0)
        scoremat[i] = mini
        arrow[i] = minpos
    return scoremat, arrow, mini


@cython.boundscheck(False)
@cython.wraparound(False)
def FastNW_(np.ndarray[DTYPEF_t, ndim=2] subvals, s1, s2, float gap=0.9):
    """Perform dynamic programing alignment of EC numbers sequences. Creates
    the scoring and arrow matrices using subvals matrix and diagonal arrays.

    Returns the dynamic programing matrix and the arrow matrix used to trace
    back the alignment.

    Arguments:
    - `subvals`: matrix of precomputed values for the alignment of each pair
                of characters
    - `s1`: list, enzymatic step sequence 1
    - `s2`: list, enzymatic step sequence 2
    - `gap`: gap penalties, default = 1
    """
    cdef  int l1, l2, li, i, j
    l1, l2 = len(s1), len(s2)
    # Create the score and arrow matrices
    cdef np.ndarray[DTYPEF_t, ndim= 2] scoremat = np.zeros((l1 + 1, l2 + 1), DTYPEF)
    cdef np.ndarray[DTYPE_t, ndim= 2] arrow = np.zeros((l1 + 1, l2 + 1), DTYPE)
    # Create first row and first column with gaps
    # scoremat[0] = np.arange(l2 + 1) * gap
    # scoremat[:, 0] = np.arange(l1 + 1) * gap
    # arrow[0] = np.ones(l2 + 1)
    for i in range(l2 + 1):
        scoremat[0, i] = i * gap
        arrow[0, i] = 1
    for i in range(l1 + 1):
        scoremat[i, 0] = i * gap
    # Compute diagonal list
    dlist = diagonal_list(l1, l2)
    # fill the matrix
    for pair in dlist:
        li = len(pair[0])
        f = np.zeros((3, li))  # results of the tree posibles xhoices
# for each  value in the diagonal
        x = pair[0]
        y = pair[1]
        f[0] = scoremat[x - 1, y] + gap
        f[1] = scoremat[x, y - 1] + gap
        f[2] = scoremat[x - 1, y - 1] + subvals[pair]
        f -= 0.001 * np.sign(f) * np.random.ranf(f.shape)  # for randomly
# select from a tie
        mini = f.min(0)
        minpos = f.argmin(0)
        scoremat[pair] = mini
        arrow[pair] = minpos
    return scoremat, arrow, mini


cdef minarg(float f[3]):
    cdef int argmin = 0
    cdef float mini = f[0]
    for i in range(1, 3):
        if f[i] < mini:
            mini = f[i]
            argmin = i
    return mini, argmin


# best performance
@cython.boundscheck(False)
@cython.wraparound(False)
def FastNW__(np.ndarray[DTYPEF_t, ndim=2] mat, dict alp, list s1, list s2,
             float gap=0.9):
    """Perform dynamic programing alignment of EC numbers sequences. Creates
    the scoring and arrow matrices using subvals matrix and diagonal arrays.

    Returns the dynamic programing matrix and the arrow matrix used to trace
    back the alignment.

    Arguments:
    - `subvals`: matrix of precomputed values for the alignment of each pair
                of characters
    - `s1`: list, enzymatic step sequence 1
    - `s2`: list, enzymatic step sequence 2
    - `gap`: gap penalties, default = 1
    """
    cdef:
        int l1, l2, li, i, j, k, eci, ecj, minpos
        float mini, randum, randfactor
        float f[3]
    mini = 0
    l1, l2 = len(s1), len(s2)
    # Create the score and arrow matrices
    cdef np.ndarray[DTYPEF_t, ndim= 2] scoremat = np.zeros((l1 + 1, l2 + 1), DTYPEF)
    cdef np.ndarray[DTYPE_t, ndim= 2] arrow = np.zeros((l1 + 1, l2 + 1), DTYPE)
    # Create first row and first column with gaps
    for i in range(l2 + 1):
        scoremat[0, i] = i * gap
        arrow[0, i] = 1
    for i in range(l1 + 1):
        scoremat[i, 0] = i * gap
    for i in range(1, l1 + 1):
        for j in range(1, l2 + 1):
            eci = alp[s1[i - 1]]
            ecj = alp[s2[j - 1]]
            sco = mat[eci, ecj]
            f[0] = scoremat[i - 1, j] + gap
            f[1] = scoremat[i, j - 1] + gap
            f[2] = scoremat[i - 1, j - 1] + sco
            # f -= 0.001 * np.sign(f) * np.random.ranf(3)  # for randomly
            # mini = f.min()
            # minpos = f.argmin()
            for k in range(3):
                randfactor = copysignf(0.001, f[k])
                randnum = rand() / float(INT_MAX)
                f[k] -= randfactor * randnum
            mini, minpos = minarg(f)
            scoremat[i, j] = mini
            arrow[i, j] = minpos
    return scoremat, arrow, mini


@cython.boundscheck(False)
@cython.wraparound(False)
def scoring(np.ndarray[DTYPEF_t, ndim=2] mat, dict alp, seq1, seq2,
            float gap=1, float fhomo=0.95, float fpengap=0.05,
            bool local=False):
    """ This function evaluates an enzymatic step sequences pair alignment.
    The sequences need to be aligned. This function tries to recall the
    fitness function in Ortegon-Cano, 2011.

    Arguments:
    - `mat`: score matrix
    - `alp`: alphabet, index of score matrix, dic
    - `s1`: enzimatic step sequence 1, list
    - `s2`: enzimatic step sequence 2, list
    - `gap`: gap penalization
    - `fhomo`: Weight factor to multiplicate the homogeneity fraction of score
    - `fpengap`: Weight factor to multiplicate the gap penalization
    - `local`: If True, evaluates the score in a local fassion, rather than
            global, i.e.  eliminates from the alignment the initial and
            final gaps of the shortest sequence and the respective positions
            in the largest sequence.
    """
    cdef unsigned int l1, l2, length, i, ii, fi, x, y, g_lenght
    cdef float homo, gaps, pengaps, score
    l1 = len(seq1)
    l2 = len(seq2)
    assert l1 == l2, 'The aligned sequences must be of the same length'
    # -- Mean column homogeneity evaluation
    if local:
        length = l1
        ii = 0
        for i in range(length):
            if seq1[i] != '-.-.-' and seq2[i] != '-.-.-':
                ii = i    # initial index
                break
        fi = length
        for i in range(length)[::-1]:
            if seq1[i] != '-.-.-' and seq2[i] != '-.-.-':
                break
            fi = i      # final index
        # trimming the sequences
        seq1 = seq1[ii:fi]
        seq2 = seq2[ii:fi]

    length = l1  # Alignment length
    homo_col = []
    for i in range(length):
        if seq1[i] == '-.-.-' or seq2[i] == '-.-.-':  # and s1[i] != s2[1]
            homo_col.append(gap)
        elif seq1[i] == '...' or seq2[i] == '...':
            pass
        elif seq1[i] == '' or seq2[i] == '':
            pass
        else:
            x = alp[seq1[i]]
            y = alp[seq2[i]]
            homo_col.append(mat[x, y])  # pondered by 0.6 as in

    homo = np.mean(homo_col)
    # -- Column increment penalization > here isnt necesary,isn't it?
    # ls1 = len([s for s in s1 if s != '-.-.-'])
    # ls2 = len([s for s in s2 if s != '-.-.-'])
    # incpen = (length - max(ls1, ls2))/10.
    # -- Gap penalization
    s1 = ':'.join(seq1)
    s2 = ':'.join(seq2)
    s1 = s1.strip('.-:')  # remove initial and final gaps
    s2 = s2.strip('.-:')

    # gap = re.compile(r"(:-.-.-)+")  # gap regex
    gaps = 0  # gap penalty value
    for s in (s1, s2):
        sgaps = []  # list of gap blocks
        for g in gapre.finditer(s):
            g = g.group()
            g_lenght = g.count("-.-.-")
            sgaps.append(g_lenght)
        if sgaps == []:
            gaps += 0
        else:
            gaps += float(len(sgaps)) / sum(sgaps)  # num of gap blocks /
# num of independent blocks per sequence
    pengap = gaps / 2  # mean , between number of sequences
    ###################################################
    # IMPORTANT ---                                   #
    # gap "penalization" argumet is not used!!!! the  #
    ###################################################
    # Final score calculation
    score = (homo * fhomo) + (pengap * fpengap)  # + incpen
    return score


def gappen(list seq1, list seq2):
    """Returns the number
    Keyword Arguments:
    list seq1 -- 
    list seq2 -- 
    """
    cdef:
        int l1 = len(seq1)
        float gappen = 0
        str gap = '-.-.-'
        int i
        float sgaps = 0
        float ngaps = 0
        float sgaps2 = 0
        float ngaps2 = 0
        float fgaps = 0
        float fgaps2 = 0
        float fblock1 = 0
        float fblock2 = 0
        str ec1, ec2
        bool begin1 = True
        bool begin2 = True
        bool block1 = True
        bool block2 = True

    # four cases depending in wich sequence has gaps
    if gap not in seq1 and gap not in seq2:
        return gappen           # no gaps, return 0
    elif gap in seq1 and gap not in seq2:
        for i in range(l1):      # Gaps only in seq1
            ec1 = seq1[i]
            if begin1 == True:
                if ec1 == gap:
                    continue
                else:
                    begin1 = False
                    continue
            else:
                if ec1 == gap:
                    sgaps += 1
                    fblock1 += 1
                    if block1 == True:
                        ngaps += 1
                        block1 = False
                else:
                    block1 = True
                    fblock1 = 0
        if block1 == False:
            sgaps -= fblock1
            ngaps -= 1
        gappen = ngaps / sgaps / 2
        return gappen
    elif gap not in seq1 and gap in seq2:
        for i in range(l1):      # Gaps only in seq2
            ec2 = seq2[i]
            if begin2 == True:
                if ec2 == gap:
                    continue
                else:
                    begin2 = False
                    continue
            else:
                if ec2 == gap:
                    sgaps += 1
                    fblock2 += 1
                    if block2 == True:
                        ngaps += 1
                        block2 = False
                else:
                    block2 = True
                    fblock2 = 0
        if block2 == False:
            sgaps -= fblock2
            ngaps -= 1
        gappen = ngaps / sgaps / 2
    else:
        for i in range(l1):
            ec1 = seq1[i]
            ec2 = seq2[i]
            if begin1 == True:
                if ec1 == gap:
                    pass
                else:
                    begin1 = False
                    pass
            else:
                if ec1 == gap:
                    sgaps += 1
                    fblock1 += 1
                    if block1 == True:
                        ngaps += 1
                        block1 = False
                else:
                    block1 = True
                    fblock1 = 0
            ec2 = seq2[i]
            if begin2 == True:
                if ec2 == gap:
                    pass
                else:
                    begin2 = False
                    pass
            else:
                if ec2 == gap:
                    sgaps2 += 1
                    fblock2 += 1
                    if block2 == True:
                        ngaps2 += 1
                        block2 = False
                else:
                    block2 = True
                    fblock2 = 0
        if block1 == False:
            sgaps -= fblock1
            ngaps -= 1
        if block2 == False:
            sgaps2 -= fblock2
            ngaps2 -= 1
        if ngaps <= 0 and ngaps2 > 0:
            gappen = (ngaps2 / sgaps2) / 2
        elif ngaps > 0 and ngaps2 <= 0:
            gappen = (ngaps / sgaps) / 2
        else:
            gappen = ((ngaps / sgaps) + (ngaps2 / sgaps2)) / 2
    return gappen

# better performance


@cython.boundscheck(False)
@cython.wraparound(False)
def scoring_(np.ndarray[DTYPEF_t, ndim=2] mat, dict alp, list seq1, list seq2,
             float gap=1, float fhomo=0.95, float fpengap=0.05):
    """ This function evaluates an enzymatic step sequences pair alignment.
    The sequences need to be aligned. This function tries to recall the
    fitness function in Ortegon-Cano, 2011.

    Arguments:
    - `mat`: score matrix
    - `alp`: alphabet, index of score matrix, dic
    - `s1`: enzimatic step sequence 1, list
    - `s2`: enzimatic step sequence 2, list
    - `gap`: gap penalization
    - `fhomo`: Weight factor to multiplicate the homogeneity fraction of score
    - `fpengap`: Weight factor to multiplicate the gap penalization
    """
    cdef unsigned int l1, l2, length, i, ii, fi, x, y, g_lenght, alisum
    cdef float gaps, pengaps, score, sgaps, ngaps
    cdef str ec1, ec2, s1, s2, s
    cdef str gapstr = '-.-.-'
    cdef str dotstr = '...'
    # define float for homogeneity mean. alisum = n aligned
    cdef float homo = 0.0
    l1 = len(seq1)
    l2 = len(seq2)
    assert l1 == l2, 'The aligned sequences must be of the same length'
    # -- Mean column homogeneity evaluation
    length = l1  # Alignment length
    # homo_col = []
    for i in range(length):
        ec1 = seq1[i]
        ec2 = seq2[i]
        if ec1 == gapstr or ec2 == gapstr:  # and s1[i] != s2[1]
            # homo_col.append(gap)
            homo += gap
        elif ec1 == dotstr or ec2 == dotstr:
            pass
        elif ec1 == '' or ec2 == '':
            pass
        else:
            x = alp[ec1]
            y = alp[ec2]
            homo += mat[x, y]  # pondered by 0.6 as in

    homo = homo / length
    pengap = gappen(seq1, seq2)
    ###################################################
    # IMPORTANT ---                                   #
    # gap "penalization" argumet is not used!!!! the  #
    ###################################################
    # Final score calculation
    score = (homo * fhomo) + (pengap * fpengap)  # + incpen
    return score


def backtrace(np.ndarray[DTYPE_t, ndim=2] arrow, list s1, list s2):
    """Reads the arrow matrix and return the aligned sequences in list form.
    Global, N-W

    Arguments:
    - `arrow`: arrow matrix generated by dynamic programing
    - `s1`: enzymatic step sequence, list
    - `s2`: enzymatic step sequence, list
    """
    #v, h = arrow.shape
    cdef:
        int v = len(s1)
        int h = len(s2)
        bool ok = True
        list st1 = []
        list st2 = []
    # Transform the sequences strings into EC numbers lists
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
    st1 = st1[::-1]
    st2 = st2[::-1]
    return st1, st2


def NW(np.ndarray[DTYPEF_t, ndim=2] mat, dict ecs, s1, s2, float gap=0.9,
       bool local=False, bool localize=False, bool nws=False):
    """
    NW alignment function. Creates a pairwise alignment using a Needelman-Wunsh
    algorithm.

    Arguments:
    - `mat`: EC number (3 levels) substitution matrix
    - `ecs`: List of ec numbers that represent the labels of the matrix
    - `s1`: EC numbers sequence 1, list
    - `s2`: EC numbers sequence 2, list
    - `gap`: gap penalty
    - `local`: if local = True, then the score of the alignment is calculated
               localy, i.e. only en the part of the alignment covered by the
               shortest sequence
    - `localize`: if localize = True, then, the function returns only the
               fragment of the alignment covered by the shortest sequence
    - `nws`: if true, the function returns the scores as is returned by the
               NW algorithm
    """
    cdef float score
    submat = FastSubValues(mat, ecs, s1, s2)  # create sub score matrix
    sco, arr, mini = FastNW_(submat, s1, s2, gap=gap)  # create score and
# arrow matrices
    s1, s2 = backtrace(arr, s1, s2)  # backtrace alingment
    # scoring : the scoring is doing using the entropy schema to
    # normalize the score in a range 0 to 1 and to generate similar
    # results to the genetic algorithm
    if nws:                    # NW score
        if localize:
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
    score = scoring(mat, ecs, s1, s2, gap=scoring_gap, local=local,
                    fhomo=0.95, fpengap=0.05)
    if localize:
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


def NW_(np.ndarray[DTYPEF_t, ndim=2] mat, dict ecs, list seq1, list seq2,
        float gap=0.9, bool local=False, bool localize=False,
        bool nws=False):
    """
    NW alignment function. Creates a pairwise alignment using a Needelman-Wunsh
    algorithm.

    Arguments:
    - `mat`: EC number (3 levels) substitution matrix
    - `ecs`: List of ec numbers that represent the labels of the matrix
    - `s1`: EC numbers sequence 1, list
    - `s2`: EC numbers sequence 2, list
    - `gap`: gap penalty for NW algorithm
    - `localize`: if localize = True, then, the function returns only the
               fragment of the alignment covered by the shortest sequence
               and the score is calculated accordingly
    - `nws`: if true, the function returns the scores as is returned by the
               NW algorithm
    """
    cdef:
        float score, mini, scoring_gap  
        str aseq1, aseq2
        int length, i, ii, fi
        char binit      # localize initial index boolena
        char bblock     # block mar boolean
    sco, arr, mini = FastNW__(mat, ecs, seq1, seq2, gap=gap)  # create score and
# arrow matrices
    seq1, seq2 = backtrace(arr, seq1, seq2)  # backtrace alingment
    # scoring : the scoring is doing using the entropy schema to
    # normalize the score in a range 0 to 1 and to generate similar
    # results to the genetic algorithm
    ### ---
    # if localize = localized alignment
    if localize:
        binit = 1
        bblock = 0
        length = len(seq1)
        ii = 0
        fi = 0
        for i in range(length):
            if seq1[i] != '-.-.-' and seq2[i] != '-.-.-' and binit:
                ii = i
                binit = 0
            elif (seq1[i] == '-.-.-' or seq2[i] == '-.-.-') and bblock:
                fi = i
                bblock = 0
            elif seq1[i] != '-.-.-' and seq2[i] != '-.-.-':
                bblock = 1
        if bblock:
            fi = length
        seq1 = seq1[ii: fi]
        seq2 = seq2[ii: fi]
    if nws:
        score = mini
    else:
        scoring_gap = 1
        score = scoring_(mat, ecs, seq1, seq2, gap=scoring_gap, fhomo=0.95,
                         fpengap=0.05)
    aseq1 = ':'.join(seq1)
    aseq2 = ':'.join(seq2)
    return aseq1, aseq2, score


        
    # if nws:                    # NW score
    #     if localize:
    #         length = len(s1)
    #         ii = 0
    #         for i in range(length):
    #             if s1[i] != '-.-.-' and s2[i] != '-.-.-':
    #                 ii = i    # initial index
    #                 break
    #         fi = length
    #         for i in range(length)[::-1]:
    #             if s1[i] != '-.-.-' and s2[i] != '-.-.-':
    #                 break
    #             fi = i      # final index
    #         s1 = s1[ii:fi]
    #         s2 = s2[ii:fi]
    #         s1 = ':'.join(s1)
    #         s2 = ':'.join(s2)
    #     return s1, s2, mini
    # scoring_gap = 1
    # score = scoring_(mat, ecs, s1, s2, gap=scoring_gap, fhomo=0.95,
    #                  fpengap=0.05)
    # if localize:
    #     length = len(s1)
    #     ii = 0
    #     for i in range(length):
    #         if s1[i] != '-.-.-' and s2[i] != '-.-.-':
    #             ii = i    # initial index
    #             break
    #     fi = length
    #     for i in range(length)[::-1]:
    #         if s1[i] != '-.-.-' and s2[i] != '-.-.-':
    #             break
    #         fi = i      # final index

    #     s1 = s1[ii:fi]
    #     s2 = s2[ii:fi]
    # aseq1 = ':'.join(s1)
    # aseq2 = ':'.join(s2)
    # return aseq1, aseq2, score

