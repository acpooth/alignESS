#
# Lines to work with alignESS.py scripts
#

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
lecs = alignESS.lecs
hmat = alignESS.hmat
db = s3.connect('nr.db')

x = db.execute("SELECT ec3 from nrseqs")
seqs_s = [i[0] for i in x]
seqs = [i.split(':') for i in seqs_s]

db.close()

#seqs = [i for i in seqs_s if len(i) > 1]
subseqs = seqs[:900]
