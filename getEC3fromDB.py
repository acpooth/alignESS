#
#
# 

""" Script for create a list form EC numbers from a enzimatic step DB
"""

import sqlite3 as s3
from numpy import sort


def EClist(db, ec='ec3', table='seqs', where='', nr=True):
    """Creates a list of ec numbers form a enzimatic step sequence DB
    
    Arguments:
    - `db`: a sqlite3 db connection
    - `ec`: level of clasifications of ec number list. ec3 or ec4. Default ec3 = ec numbers 3th level of clasification
    - `table`: table used to generate the list
    - `where`: a sqlite where statement
    - `nr`: if nr=True, returns a no redundant list of ECs, else, returns a list with all the times an EC appears in the db 
    """
    qwery = db.execute("""select {} from {} {}""".format(ec, table, where)).fetchall()
    qwery = [q[0] for q in qwery]
    ecs = [ec for seq in qwery for ec in seq.split(':')]
    if nr:
        ecs = sort(list(set(ecs)))
    else:
        ecs = sort(ecs)
    return tuple(ecs)


