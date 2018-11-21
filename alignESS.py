#!/usr/bin/python3
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
"""Program to align Enzimatic Step Sequences (ESS) using Dynamic Programming (DP)
and Genetic Algorithms (GA)"""

import os
import sys
import argparse
import numpy as np
import sqlite3 as s3
import pyximport
pyximport.install(setup_args={"include_dirs": np.get_include()},
                  reload_support=True)
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
    sep   -- String that separates fields (default TAB)
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
        esslist = []
        indexlist = []
        for es in esss:
            if es == '':
                continue
            es = es.split(sep)
            esslist.append(es[1].split(':'))
            indexlist.append(es[0])
        return esslist, indexlist


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


def alignESS(seq1, seq2, gap=0.9, localize=False):
    """Wrapper for nwx.NW function for ESS alignment

    Keyword Arguments:
    seq1 -- ESS 1. str, colon separated EC numbers.
    seq2 -- ESS 1. str, colon separated EC numbers.
    gap  -- Gap penalization.
    localize -- if localize = True, the funtion returns only
            the fragment of the alignment covered by the
            shortest sequence and the score is calculated
            accordingly
    """
    seq1 = seq1.split(':')
    seq2 = seq2.split(':')
    seq1, seq2, score = nwx.NW(
        hmat, decs, seq1, seq2, gap=gap, localize=localize)
    return seq1, seq2, score

##########
# Parser #
##########


def arg_parser(get_parser=False):
    """If parser == True,  the program retunr parser obhect"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-m', '--matrix', default='h_ent_mat.npz',
                        help='''EC number similarity matrix
                        (TO DO)[h_ent_mat.npz]''')
    parser.add_argument('--gap', type=float, default=0.9,
                        help='Gap penalization (from 0 to 1) [0.9]')
    # TO DO parser -c, --create-matrix,
    # 'create matrix given a list of ec numbers'
    ##############
    # Subparsers #
    ##############
    subparsers = parser.add_subparsers(help='ESSs alignment subcommands',
                                       dest='command')
    # Pair alignment
    pairp = subparsers.add_parser('pair',
                                  help='''ESS command line pairwise
                                  comparisson using DP''')
    pairp.add_argument('ess1', type=str,
                       help='''ESS (3 levels EC numbers). Colon separated.
                        (1.2.3:3.5.-:...:9.9.9)''')
    pairp.add_argument('ess2', type=str,
                       help='''ESS (3 levels EC numbers). Colon separated.
                        (1.2.3:3.5.-:...:9.9.9)''')
    pairp.add_argument('-l', '--localize', action='store_true', default=False,
                       help='''The alignment is trimmed to the coverage of the
                       shortest ESS and the score is then calculated to the
                       trimmed alignment. This method allows to find
                       'local-like' alignments between ESS of different
                       size''')
    # DB alignment
    dbp = subparsers.add_parser('dbalign',
                                help='ESSs database(s) alignment using DP')
    dbp.add_argument('essdb1', type=str,
                     help='''ESSs database 1. Sqlite3 file with nrseqs table
                     or text file with one ESS in each line. If the essdb2
                     argument is not specified, the program performs the
                     all-vs-all alignment in essdb1. This argument
                     also can be a single ESS, in this case the -db2 argument
                     is necessary''')
    dbp.add_argument('-db2', '--essdb2', type=str,
                     help='''ESSs databse 2. Sqlite3 file with nrseqs table
                     or text file with one ESS in each line''')
    dbp.add_argument('-o', '--outfile', type=str,  default='output.txt',
                     help="""Outfile name to report scores. By default the
                     file only contains the id of the ESSs and the score
                     of the alignment. If argument '-align' is set,
                     then the file contains the alignmed ESSs""")
    dbp.add_argument('-t', '--threshold', type=float, default=0.3,
                     help='''Threshold score to filter results in the
                     range 0-1 [0.3]. If the threshold is high (>0.6) and
                     the databases are large, results may saturate the
                     RAM memmory, beware!''')
    dbp.add_argument('-nproc', type=int, default=2,
                     help='''Number of processes to execute analysis [2].
                     It can be created more processes than cores in the
                     the processor, so the speedup of the analysis
                     depends on the number of cores available''')
    dbp.add_argument('-l', '--localize', action='store_true', default=False,
                     help='''The alignment is trimmed to the coverage of the
                     shortest ESS and the score is then calculated to the
                     trimmed alignment. This method allows to find
                     'local-like' alignments between ESS of different
                     size''')
    dbp.add_argument('-align', action='store_true', default=False,
                     help='''If set, the outputfile contains the alignment
                     of each ESS pair bellow the threshold. Beware, if the
                     databases are large and the threshold high,
                     the file may be huge or the RAM memmory colapse.''')
    # Multiple Alignment
    multip = subparsers.add_parser('multi',
                                   help='ESSs multiple alignment using GA')
    multip.add_argument('multifile', type=str, metavar='FILENAME',
                        help='''ESSs file. Each line must contain an ESS name
                        and the ESS separeated by a TAB.
                        Accepts commentaries with '#' ''')
    multip.add_argument('-o', '--multiout', type=str, metavar='OUTPUTFILE',
                        default='multiout.txt',
                        help='''Multiple alignment outputfile''')
    multip.add_argument('-p', '--pcomp', type=str, metavar='FILENAME',
                        help='''If spificied, stores the pairwise comparisson
                        of the ESSs in the ESSs file''')
    # -
    if get_parser:
        return parser
    args = parser.parse_args()
    return args

################
# Main helpers #
################


def _ess_type(input_str):
    """Returns the type of ESSs recived by command line
    ['ess', 'sqdb', 'textdb']

    Keyword Arguments:
    input_str -- string from command line
    """
    if os.path.exists(input_str):  # is it a file?
        if '.db' in input_str:
            answer = 'sqlite'
        elif '.txt' in input_str:
            with open(input_str) as inf:
                line = inf.readline()
            line = line.split('\t')
            if len(line) == 1:
                answer = 'text_noind'
            elif len(line) >= 2:
                answer = 'text_ind'
        else:
            answer = 'file'
    else:                       # is it a ESS?
        # asuming the count of : and .
        # for each : you need 2(.) + 2
        colons = input_str.count(':')
        points = input_str.count('.')
        if (colons >= 1
                and points == colons * 2 + 2):
            answer = 'ess'
        else:
            print('File not found or ESS in wrong format')
            print('Exiting')
            raise FileNotFoundError('''"{}" file not found or wrong ESS
format.'''.format(input_str))
    return answer


def _loaddb(fname):
    """Load a database according to type

    Keyword Arguments:
    fname -- filename
    """
    _type = _ess_type(fname)
    if _type == 'sqlite':
        return load_sqlite(fname)
    elif _type == 'text_ind':
        return load_txt(fname, index=True, sep='\t')
    elif _type == 'text_noind':
        return load_txt(fname), None
    elif _type == 'file':
        return load_txt(fname)
    else:
        print('Something is really wierd')


def _load_multi(fname):
    """Loads a multiple alignment filename.
    2 columns tabular format:
    Col 1 = ESS name
    Col 2 = ESS 

    accept commentaries with #

    Keyword Arguments:
    fname -- filename
    """
    names = []
    esss = []
    with open(fname) as inf:
        for line in inf:
            line = line.strip()
            if not line:
                continue
            if line[0] == '#' or line[0] == '>':
                continue
            name, ess = line.split('\t')
            ess = ess.split(':')
            names.append(name)
            esss.append(ess)
    return esss, names


##################
# Main functions #
##################


def main_pair(args):
    """Main function for pair subcommand
    """
    from sys import stdout
    seq1, seq2, score = alignESS(args.ess1, args.ess2, gap=args.gap,
                                 localize=args.localize)
    stdout.write('ess1:\t{}\n'.format(seq1))
    stdout.write('ess2:\t{}\n'.format(seq2))
    stdout.write('score = {}\n'.format(score))


def main_db(args):
    """Main function for dbalign subcomand
    """
    print('------ Opening databases:')
    typ1 = _ess_type(args.essdb1)
    print(args.essdb1, '(', typ1, ')')
    if args.essdb2:
        typ2 = _ess_type(args.essdb2)
        assert typ2 != 'ess', '-db2 must be a database, not an ESS'
        print(args.essdb2, '(', typ2, ')')
    if args.align is True:
        oscore = False
    else:
        oscore = True
    # Checking types
    if not args.essdb2:
        if typ1 == 'ess':
            print('#' * 20)
            print('Warning!!! ---------')
            print('When using an ESS as input, you need to spicify ' +
                  'an ESS database in -db2 argument.')
            print('Exit program!!! :D try again!!!')
            exit()
        db1, indices = _loaddb(args.essdb1)
        print('------ Aligning the database with itself (all vs all)...')
        print('Number of processes: {}'.format(args.nproc))
        scores = nwx.alldb_comp(db1, hmat, decs, thres=args.threshold,
                                nproc=args.nproc,
                                localize=args.localize,
                                oscore=oscore
                                )
        print('------ Storing data:', args.outfile)
        nwx.store_dict(args.outfile, scores, indices=indices, indices2=indices)
        return
    if typ1 == 'ess':
        db2, ind2 = _loaddb(args.essdb2)
        print('------ Aligning ESS vs DB ------')
        print('Number of processes: 1 (fixed)')
        ess = args.essdb1.split(':')
        scores = nwx.seq_vs_db(ess, db2, hmat, decs, thres=args.threshold,
                               localize=args.localize,
                               oscore=oscore
                               )
        print('------ Storing data:', args.outfile)
        print('Number of hits: ', len(scores))
        if not scores:
            print('WARNING: 0 hits, file {} not created'.format(args.outfile))
        else:
            nwx.store_dict(args.outfile, scores, indices=ind2)
    else:
        # if two databases were speficied
        db1, ind1 = _loaddb(args.essdb1)
        db2, ind2 = _loaddb(args.essdb2)
        print('------ Aligning both databases ------')
        print('Number of processes: {}'.format(args.nproc))
        scores = nwx.db_vs_db(db1, db2, hmat, decs, thres=args.threshold,
                              nproc=args.nproc,
                              localize=args.localize,
                              oscore=oscore
                              )
        print('------ Storing data:', args.outfile)
        nwx.store_dict(args.outfile, scores, indices=ind1, indices2=ind2)


def main_multi(args):
    """Main function for multiple alignment
    """
    # pcomp
    import tempfile
    from subprocess import run, PIPE
    print('------ Loading file')
    print(args.multifile)
    esss, names = _load_multi(args.multifile)
    print('------ Making pairwise comparissons')
    pairs = nwx.alldb_comp(esss, hmat, decs)  # pair align
    if args.pcomp:
        nwx.store_dict('args.pcomp', pairs)
    # create temp files
    tempsco = tempfile.NamedTemporaryFile(mode='w+t', dir='.')
    templist = tempfile.NamedTemporaryFile(mode='w+t', dir='.')
    # fill files
    for i, nam in enumerate(names):
        ess = ':'.join(esss[i])
        line = "{}\t{}\t{}\n".format(i, nam, ess)
        templist.write(line)
    for i, dic_ in pairs.items():
        for j, sco in dic_.items():
            line = "{}\t{}\t{}\n".format(i, j, sco)
            tempsco.write(line)
    # rewind files
    tempsco.seek(0)
    templist.seek(0)
    # Run multiple alignment
    binf = 'bin/AlineaMultiple'
    binpath = os.path.join(exedir, binf)
    npob, ngen, cross, mut = '100', '200', '0.7', '0.1'
    pengap, homo, peninc = '0.05', '0.9', '0.05'
    cmd = [binpath, templist.name, tempsco.name,
           npob, ngen, cross, mut, pengap, homo, peninc, '1']
    print('------ Building multiple alignment')
    print("""- Genetic algorithm parameters:
Population: {}
Max generations: {}
Crossover prob: {}
Mutation prob: {}
- Objetive function parameters (must sum 1):
Gap penalization: {}
Homogeneity: {}
Column increment penalization:{}
""".format(npob, ngen, cross, mut, pengap, homo, peninc))
    malign = run(cmd, stdout=PIPE)
    malign = malign.stdout.decode('utf8')
    # End temp files
    tempsco.close()
    templist.close()
    print('------ Creating file')
    print(args.multiout)
    with open(args.multiout, 'w') as outf:
        for line in malign.split('\n'):
            if line == '':
                continue
            elif line[0] == 'F':
                newl = '# ' + line + '\n'
            else:
                line = line.split('\t')
                newl = '{}\t{}\n'.format(line[1], line[2])
            outf.write(newl)


def main():
    args = arg_parser()
    if args.command == 'pair':
        main_pair(args)
    elif args.command == 'dbalign':
        main_db(args)
    elif args.command == 'multi':
        main_multi(args)
    elif args.command is None:
        parser = arg_parser(True)
        parser.print_help()
    print('>>> Done!!! <<<')
    print(':D, see you soon.')


if __name__ == '__main__':
    main()
