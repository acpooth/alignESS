from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("dlist.pyx")
)

%load_ext autoreload
%load_ext Cython
%autoreload 2
import nw_ec_align as nw
lecs = nw.lecs
decs = nw.decs

import pyximport
pyximport.install()
import dlist
import nw_ec_alignx as nwx


# subval = nw.__FastSubValues(nw.hmat, decs, nw.ss1, nw.ss2)
# align = dlist.FastNW(subval, nw.ss1, nw.ss2)
# arrow = align[1]
# arrowint = arrow.astype(np.int)
# back_ = dlist.backtrace(arrowint, nw.ss1, nw.ss2)
# back = nw.backtrace(arrow, nw.s1, nw.s2)

%timeit - n 100 nw.NW(nw.hmat, lecs, nw.s1, nw.s2)
%timeit - n 100 dlist.NW_(nw.hmat, decs, nw.ss1, nw.ss2)
%timeit - n 100 nwx.NW(nw.hmat, decs, nw.ss1, nw.ss2)

ess = [nw.ss1, nw.ss2, nw.ss3, nw.ss1_]
print('\tOld\t|\tNew')
print('Sco\tLoc\t|Sco\tLoc')
print('---------------------')

line = '{:.3}\t{:.3}\t|{:.3}\t{:.3}'
for i in ess:
    for j in ess:
        sold = nw.NW(nw.hmat, lecs, ':'.join(i), ':'.join(j))[2]
        snew = nwx.NW(nw.hmat, decs, i, j)[2]
        soldl = nw.NW(nw.hmat, lecs, ':'.join(i), ':'.join(j),
                      local=True, localize=True)[2]
        snewl = nwx.NW(nw.hmat, decs, i, j, localize=True)[2]
        print(line.format(sold, soldl, snew, snewl))


###############################
# Error divition by 0         #
###############################
import sqlite3 as s3
import alignESS

db = s3.connect('nr.db')

x = db.execute("SELECT ec3 from nrseqs")
seqs = [i[0].split(':') for i in x]
seqs = [i for i in seqs if len(i) > 1]

mat = alignESS.hmat
decs = alignESS.decs
as1 = ['1.1.1', '1.1.1', '1.1.1', '4.1.3', '4.2.1']
as2 = ['1.1.1',  '1.1.1', '-.-.-', '-.-.-', '-.-.-']
alignESS.nwx.scoring(mat, decs, as1, as2)


ss1 = '1.1.1:1.1.1:1.1.1:4.1.3:4.2.1:-.-.-:-.-.-'.split(':')
ss2 = '-.-.-:1.1.1:2.6.1:4.1.2:4.2.3:2.7.1:1.1.1'.split(':')

ss1 = ['1.1.1', '1.1.1', '1.1.1', '-.-.-', '4.1.3', '4.2.1', '-.-.-', '-.-.-']
ss2 = ['-.-.-', '-.-.-', '1.1.1', '2.6.1', '4.1.2', '4.2.3', '2.7.1', '1.1.1']
alignESS.nwx.scoring(mat, decs, as1, as2)


def pairali(s1, s2): return alignESS.nwx.NW(
    alignESS.hmat, alignESS.decs, s1, s2)


error = 0
for i in seqs:
    try:
        s1, s2, sco = pairali(i, seqs[15])
        if sco <= 0.3:
            print(sco)
    except:
        error += 1
print(error)

# number of error reduced by adding EC 9.9.9 to the sim matrix
# number of errors is reduced by modifing gappen function for 0 division

erseq = []
err = []
for i in seqs:
    try:
        s1, s2, sco = pairali(i, seqs[15])
        if sco <= 0.3:
            print(sco)
    except (ZeroDivisionError) as e:
        erseq.append(i)
        err.append(e)


arr, mini = alignESS.nwx._FastNW(mat, decs, seqs[15], erseq[0])
as1, as2 = alignESS.nwx._backtrace(arr, seqs[15], erseq[0])
alignESS.nwx.scoring(mat, decs, as1, as2)

for seq in erseq:
    ar_, min_ = alignESS.nwx._FastNW(mat, decs, seqs[15], seq)
    _as1, _as2 = alignESS.nwx._backtrace(ar_, seqs[15], seq)
    try:
        alignESS.nwx.scoring(mat, decs, _as1, _as2)
    except:
        print(':'.join(_as1))
        print(':'.join(_as2))


from random import choice


for i in range(100):
    seq1 = choice(seqs)
    print('Working with: ', ':'.join(seq1))
    for seq2 in seqs[:8000]:
        pairali(seq1, seq2)


########################
# Multipthreading test #
########################

import sqlite3 as s3
import alignESS

decs = alignESS.decs
mat = alignESS.hmat

db = s3.connect('nr.db')

x = db.execute("SELECT ec3 from nrseqs")
seqs = [i[0].split(':') for i in x]
seqs = [i for i in seqs if len(i) > 1]

subseqs = seqs[:900]

# Multiprocessing


def aliall(ind, seqs, mat, decs):
    """Align the sequences in the specified index (ind) versus
all the rest of sequences
    Keyword Arguments:
    ind  --
    seqs --
    """
    seq1 = seqs[ind]
    result = []
    for i in range(ind + 1, len(seqs)):
        seq2 = seqs[i]
        s1, s2, sco = alignESS.nwx.NW(mat, decs, seq1, seq2)
        result.append(sco)
    return result


from multiprocessing import Pool, Array
from functools import partial


align_func = partial(aliall, mat=mat, decs=decs, seqs=subseqs)
pool = Pool(processes=2)
indices = list(range(500))
res = pool.map(align_func, indices, chunksize=1000)

i1, i2 = np.triu_indices(total)
i1 = np.int16(i1 + 1)
i2 = np.int16(i2 + 1)
indices = np.vstack((i1, i2)).T
# scomat = np.zeros((total, total), dtype = np.float16)
# alternative scomat
shared_sco_base = Array(ctypes.c_float, total * total)
# scomat = np.ctypeslib.as_array(shared_sco_base.get_obj())
scomat_ = np.frombuffer(shared_sco_base.get_obj(), dtype=np.float32)
scomat = scomat_.reshape((total, total))
global scomat
align_func = partial(fill_mat,  nrdic=nrdic)
pool = Pool(processes=6)
pool.map(align_func, indices, chunksize=1000)
# ids1 and ids2 are not strictly necessary because data is stored in


##########################################################
# Using cython paralell.prange does not working          #
# because in order to execute the NW function is         #
# necessary to activate the GIL. This is necessary       #
# because this function contains a lot of python         #
# variables :S                                           #
#                                                        #
# This may function if the dinamic programming algorithm #
# uses only c variables.                                 #
##########################################################

pairall = alignESS.nwx.pairall
pairall_t = alignESS.nwx.pairall_t

pairall(subseqs, decs, mat)
pairall_t(subseqs, decs, mat)

%timeit - n 3 pairall(subseqs, decs, mat)
%timeit - n 3 pairall_t(subseqs, decs, mat)


####################################
# Threading test -- not working... #
# you forgot the GIL!!!            #
####################################


# import threading
# def aliall(ind, seqs, mat, decs):
#     """Align the sequences in the specified index (ind) versus
# all the rest of sequences
#     Keyword Arguments:
#     ind  --
#     seqs --
#     """
#     seq1 = seqs[ind]
#     result = []
#     # print('Thread starting')
#     for i in range(ind + 1, len(seqs)):
#         seq2 = seqs[i]
#         s1, s2, sco = alignESS.nwx.NW(mat, decs, seq1, seq2)
#         result.append(sco)
#     print('Thread ending: ', len(result))
#     # print(result[:10])
#     # print(len(result))
#     return result


# %%timeit -n 4
# for i in range(4):
#     aliall(0, seqs, mat, decs)

# %%timeit -n 4
# threads = []
# for i in range(4):
#     t =threading.Thread(target=aliall, args=(0, seqs, mat, decs, ))
#     t.start()
#     threads.append(t)

# x = threading.Thread(target=aliall, args=(0, seqs, mat, decs, ))

# x = thread.start_new_thread(aliall, (0, seqs, ))
# y = thread.start_new_thread(aliall, (0, seqs, ))
# z = thread.start_new_thread(aliall, (0, seqs, ))
# zz = thread.start_new_thread(aliall, (0, seqs, ))


# thread.start_new_thread(print_time, ('Thread-3', 8, ))
