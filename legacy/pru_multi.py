

from sys import argv
from multiprocessing import Pool, Array
from functools import partial
import time
import sqlite3 as s3
import alignESS
import nw_ec_align as nw
import pyximport
pyximport.install(reload_support=True)
import nw_ec_alignx as nwx


decs = alignESS.decs
mat = alignESS.hmat
db = s3.connect('nr.db')

x = db.execute("SELECT ec3 from nrseqs")
seqs_s = [i[0] for i in x]
seqs = [i.split(':') for i in seqs_s]

#seqs = [i for i in seqs_s if len(i) > 1]
subseqs = seqs[:900]


def aliall(ind, seqs, mat, decs):
    """Align the sequences in the specified index (ind) versus
all the rest of sequences
    Keyword Arguments:
    ind  --
    seqs --
    """
    seq1 = seqs[ind]
    result = []
    print(ind)
    for i in range(ind + 1, len(seqs)):
        seq2 = seqs[i]
        s1, s2, sco = alignESS.nwx.NW(mat, decs, seq1, seq2)
        result.append(sco)
    return result


def aliall_d(ind, seqs, mat, decs):
    """Align the sequences in the specified index (ind) versus
all the rest of sequences
    Keyword Arguments:
    ind  --
    seqs --
    """
    resdic = {ind: {}}
    seq1 = seqs[ind]
    print(ind)
    for i in range(ind + 1, len(seqs)):
        seq2 = seqs[i]
        s1, s2, sco = alignESS.nwx.NW(mat, decs, seq1, seq2)
        resdic[ind][i] = sco
    return resdic


def aliall_old(ind, seqs, mat, ecs):
    """Align the sequences in the specified index (ind) versus
all the rest of sequences
    Keyword Arguments:
    ind  --
    seqs --
    """
    seq1 = seqs[ind]
    result = []
    print(ind)
    for i in range(ind + 1, len(seqs)):
        seq2 = seqs[i]
        s1, s2, sco = nw.NW(mat, ecs, seq1, seq2)
        result.append(sco)
    return result



# -- multi
nprocs = int(argv[1])
nind = int(argv[2])
pool = Pool(processes=nprocs)
indices = list(range(nind))


print('Nproces: ', nprocs, ' --- nindices: ', nind)

# - lists
t1 = time.time()
print("Cython function - list")
align_func = partial(aliall, mat=mat, decs=decs, seqs=seqs)
resl = pool.map(align_func, indices, chunksize=nprocs)
resl = [i for L in resl for i in L]
t2 = time.time()
tlist = t2 - t1
del resl

# - dict FASTEST (simple python)
t1 = time.time()
print("Cython function - dict")
align_func = partial(aliall_d, mat=mat, decs=decs, seqs=seqs)
resd = pool.map(align_func, indices, chunksize=nprocs)
result = {}
for d in resd:
    result.update(d)
del resd
t2 = time.time()
tdic = t2 - t1

# # - old
# t1 = time.time()
# print("Pure python funct")
# align_func = partial(aliall_old, mat=nw.hmat, ecs=nw.lecs, seqs=seqs_s)
# res = pool.map(align_func, indices, chunksize=nprocs)
# t2 = time.time()
# told = t2 - t1


print('Time cython list = ', tlist)
print('Time cython dict = ', tdic)
# print('Time old = ', told)


# print(len(res))
# print(len(res[0]))
# print(len(res[-1]))

# Raw speed --- No storing data
# $ python3 pru_multi.py 4 200
# Time new =  50.09322452545166
# Time old =  1704.5306015014648


# Try dictionary MANAGER
# it does not work because the manager needs to unlock the
# object before using it....

# if there is no manager, the modifications in de dictionary
# are not stored
# from sys import argv
# from multiprocessing import Pool, Manager
# from functools import partial
# import time
# import sqlite3 as s3
# import alignESS
# import nw_ec_align as nw


# decs = alignESS.decs
# mat = alignESS.hmat
# db = s3.connect('nr.db')

# x = db.execute("SELECT ec3 from nrseqs")
# seqs_s = [i[0] for i in x]
# seqs = [i.split(':') for i in seqs_s]

# manager = Manager()
# resdic = manager.dict()


# def aliall(ind, seqs, mat, decs, resdic):
#     """Align the sequences in the specified index (ind) versus
# all the rest of sequences
#     Keyword Arguments:
#     ind  --
#     seqs --
#     """
#     seq1 = seqs[ind]
#     print(ind)
#     for i in range(ind + 1, len(seqs)):
#         seq2 = seqs[i]
#         s1, s2, sco = alignESS.nwx.NW(mat, decs, seq1, seq2)
#         if ind not in resdic:
#             resdic[ind] = {i: sco}
#         else:
#             resdic[ind][i] = sco
#     # return result


# nprocs = 2
# nind = 10
# pool = Pool(processes=nprocs)
# indices = list(range(nind))


# print(" Cython function")
# align_func = partial(aliall, mat=mat, decs=decs, seqs=seqs, resdic=resdic)
# res = pool.map(align_func, indices, chunksize=nprocs)


######
#  Cython 'optimiztion ' The code runs almost at the same speed if no variables declared
#####

# def alldb_shared(seqs,  mat,  decs,
#                  thres=1.0, localize=False, nproc=2):
#     """Align all the ESS in the database (seqs) and return the results in form
#     of dictionary. Thre result is equivalent to the upper part the all vs all
#     comparisson matrix.

#     BEWARE if the database is huge and the threshold is > 0.4 the procces may
#     use all the RAM

#     Keyword Arguments:
#     seqs     --
#     mat      --
#     decs     --
#     thres    -- (default 1.0)
#     localize -- (default False)
#     nproc    -- (default 2)
#     """
#     total = len(seqs)
#     create_global_arr(total)
#     indices = range(total)
#     pool = Pool(processes=nproc)
#     align_func = partial(fill_mat, seqs=seqs, mat=mat, decs=decs,
#                          thres=thres, localize=localize)
#     pool.map(align_func, indices, chunksize=nproc)
#     return scomat


# def create_global_arr(int total):
#     global scomat
#     shared_sco_base = Array(ctypes.c_float, total * total)
#     scomat_ = np.frombuffer(shared_sco_base.get_obj(), dtype=np.float32)
#     scomat = scomat_.reshape((total, total))


# @cython.boundscheck(False)
# @cython.wraparound(False)
# def fill_mat(int ind, list seqs, np.ndarray[DTYPEF_t, ndim=2] mat,
#              dict decs, bool localize, float thres):
#     cdef:
#         int j
#         list seq1, seq2
#         float sco
#         float[:, :] scomatrix = scomat
#     seq1 = seqs[ind]
#     for j in range(ind + 1, len(seqs)):
#         seq2 = seqs[j]
#         sco = NW(mat, decs, seq1, seq2, oscore=True,
#                  localize=localize)
#         if sco <= thres:
#             scomatrix[ind, j] = sco
