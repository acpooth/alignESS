#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:        nw_ec_align.px
# Purpose:     Program to align Enzimatic Step Sequences (ESS) using dynamic
#              programming
#
# @uthor:      acph - dragopoot@gmail.com
#
# Created:     Sep, 2017
# Copyright:   (c) acph 2017
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------
"""Program to align Enzimatic Step Sequences (ESS) using dynamic programming"""

import os
import argparse
import numpy as np
import sqlite3 as s3
import pyximport
pyximport.install()
import nw_ec_alignx as nwx


# load similarity matrix
matfile = 'h_ent_mat.npz'       # ecs_entropy.py generated file
# exedir = os.path.dirname(argv[0])
# exepath = os.path.join(exedir, matfile)
exedir = os.path.realpath(__file__)
exedir = os.path.dirname(exedir)
exepath = os.path.join(exedir, matfile)

hmat = np.load(exepath)
ecs = hmat['ecs']
hmat = hmat['matrix']
lecs = [str(ec) for ec in ecs]
decs = {}
for i, v in enumerate(lecs):
    decs[v] = i

# -- End load matrices

##################
# Load databases #
##################


def load_txt(fname, index=False, sep='\t'):
    """Loads a database of ESS stored in a text file. Unless index = True, each
    line correspond to one ESSs.

    Returns a list of ESS or a dictionary. Each ESS is a list of EC numbers

    Keyword Arguments:
    fname -- String, File name
    index -- If True, each line is splited according to sep (default False).
             The first element is trated as index and the second as the ESS.
             If the line contains more elements, will be ignored.
    sep   -- (default '\t')
    """
    if not index:
        with open(fname) as inf:
            esss = inf.read()
        esss = esss.split('\n')
        assert len(esss[0].split()) == 1, 'Lines must contain only the ESS'
        esss = [es.split(':') for es in esss if es != '']
        return esss
    else:
        with open(fname) as inf:
            esss = inf.read()
        esss = esss.split('\n')
        assert len(esss[0].split(sep)) >= 2, 'Lines must contain index and ESS'
        essdic = {}
        for es in esss:
            if es == '':
                continue
            es = es.split(sep)
            essdic[es[0]] = es[1].split(':')
        return essdic


def load_sqlite(dbf, table='nrseqs', length=1):
    """Loads a database of ESS sotored in a sqlite3 database.

    Returns a list of ESSs where each ESS is a list of EC numbers
    and a list of indices.

    Keyword Arguments:
    dbile  -- String, sqlite3 database filename
    table  -- String, database table name (default 'nrseqs')[seqs|nrseqs]
    length -- Int, minimum longitude of the ESS to retrive (default 1)
    """
    assert table in ['nrseqs', 'seqs'], "Table name not known [nrseqs|seqs]"
    if table == 'nrseqs':
        idname = 'nrid'
    else:
        idname = 'id'
    ids = []
    esss = []
    with s3.connect(dbf) as db:
        s3state = "SELECT {}, ec3 from {}".format(idname, table)
        cursor = db.execute(s3state)
        for entry in cursor:
            ess = entry[1].split(':')
            if len(ess) >= length:
                ids.append(entry[0])
                esss.append(ess)
    return esss, ids


######################
# Alignment wrappers #
######################


def alignESS(seq1, seq2, gap=0.9):
    """Wrapper for nwx.NW function for ESS alignment

    Keyword Arguments:
    seq1 -- ESS 1. str, colon separated EC numbers.
    seq2 -- ESS 1. str, colon separated EC numbers.
    gap  -- Gap penalization.
    """
    seq1 = seq1.split(':')
    seq2 = seq2.split(':')
    seq1, seq2, score = nwx.NW(hmat, decs, seq1, seq2, gap=gap)
    return seq1, seq2, score


def arg_parser(get_parser=False):
    """If parser == True,  the program retunr parser obhect"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-m', '--matrix', default='h_ent_mat.npz',
                        help='EC number similarity matrix (TO DO)[h_ent_mat.npz]')
    parser.add_argument('--gap', type=float, default=0.9,
                        help='Gap penalization (from 0 to 1) [0.9]')
    # parser '-c', '--create-matrix', 'create matrix given a list of ec numbers'
    ##############
    # Subparsers #
    ##############
    subparsers = parser.add_subparsers(help='ESSs alignment subcommands',
                                       dest='command')
    # Pair alignment
    pairp = subparsers.add_parser('pair',
                                  help='ESS command line pairwise comparisson')
    pairp.add_argument('ess1', type=str,
                       help='''ESS (3 levels EC numbers). Colon separated.
                        (1.2.3:3.5.-:...:9.9.9)''')
    pairp.add_argument('ess2', type=str,
                       help='''ESS (3 levels EC numbers). Colon separated.
                        (1.2.3:3.5.-:...:9.9.9)''')
    # DB alignment
    dbp = subparsers.add_parser('dbalign',
                                help='ESSs database(s) alignment')
    dbp.add_argument('essdb1', type=str,
                     help='''ESSs database 1. Sqlite3 file with nrseqs table
                     or text file with one ESS in each line. If the essdb2
                     argument is not specified, the program performs the
                     all-vs-all alignment in essdb1. This argument
                     also can be a sigle ESS''')
    dbp.add_argument('-db2', '--essdb2', type=str,
                     help='''ESSs databse 2. Sqlite3 file with nrseqs table
                     or text file with one ESS in each line''')
    dbp.add_argument('outfile', type=str,
                     help="""Outfile name to report scores. By default the
                     file only contains the id of the ESSs and the score
                     of the alignment. If argument '-align' is set, then 
                     the file contains the alignment of ESSs""")
    dbp.add_argument('-t', '--threshold', type=float, default=0.6,
                     help='''Threshold score to filter results in the
                     range 0-1 [0.3]. If the threshold is high (>0.6) and
                     the databases are large, results may saturate the
                     RAM memmory, beware!''')
    dbp.add_argument('-nproc', type=int, default=2,
                     help='''Number of processes to execute analysis [2].
                     It can be created more processes than cores in the
                     the processor, so the speedup of the analysis
                     depends on the number of cores available''')
    dbp.add_argument('-align', action='store_true',
                     help='''If set, the outputfile contains the alignment
                     of each ESS pair bellow the threshold. Beware, if the
                     databases are large, the file may be huge''')
    # -align argument
    # Multiple Alignment
    multip = subparsers.add_parser('multi',
                                   help='ESSs multiple alignment')
    # -
    if get_parser:
        return parser
    args = parser.parse_args()
    return args

##################
# Main functions #
##################


def main_pair(args):
    """Main function for pair subcommands

    """
    from sys import stdout
    seq1, seq2, score = alignESS(args.ess1, args.ess2, gap=args.gap)
    stdout.write('ess1:\t{}\n'.format(seq1))
    stdout.write('ess2:\t{}\n'.format(seq2))
    stdout.write('score = {}\n'.format(score))


def main():

    args = arg_parser()
    if args.command == 'pair':
        main_pair(args)
    elif args.command == 'dbalign':
        print(args.command)
    elif args.command == 'multi':
        print(args.command)
    elif args.command is None:
        parser = arg_parser(True)
        parser.print_help()


if __name__ == '__main__':
    main()
