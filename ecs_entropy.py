#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
# Name:     ecs_entropy.py
# Purpose:  Create an entropy based similarity (distance) matrix for EC numbers
#           alignment.
#
# @uthor:      acph - dragopoot@gmail.com
#
# Created:     Aug 2017
# Copyright:   (c) acph 2017
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#------------------------------
"""Create an entropy based similarity (distance) matrix for EC alignment.
The similarity is calculated for the first 3 levels of EC classification"""

import argparse
import numpy as np


def parser():
    """Command line argument parser"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ecfile',
                        help="""KEGG EC numbers list.
                        This file can be downloaded with this link: 
                        http://rest.kegg.jp/list/ec""")
    parser.add_argument('-s', '--store-format', default='npz',
                        choices=['npz', 'txt'], metavar='npz|txt',
                        help='''Output file format. The default (npz)
                        format is as numpy binary file that can be used by 
                        the ESS alignment program. txt option gets a tabular 
                        separated text file''')
    parser.add_argument('-p', '--plot', action='store_true',
                        help='''Plot the matrix in ent_matrix.png file''')
    args = parser.parse_args()
    return args


def load_ec_file(ecfile, ecl=3):
    """Loads the EC numbers in the ecfile into a python list

    Keyword Arguments:
    ecfile -- KEGG EC numbenrs list.
            This file can be downloaded with this link:
            http://rest.kegg.jp/list/ec
    ecl    -- Level of EC number classification. 2|[3]|4 .
            Level of EC number classification to retrive.
    """
    assert ecl in [2, 3, 4], "'ecl' must be int 2,3 or 4"
    with open(ecfile) as inf:
        ecs = []
        for line in inf:
            ec = line.split()[0]
            ec = ec.split(':')[1]
            ecs.append(ec)
    if ecl == 4:
        ecs.sort()
        return ecs
    lecs = set([])
    for ec in ecs:
        ec = ec.split('.')[:ecl]
        ec = '.'.join(ec)
        lecs.add(ec)
    lecs = list(lecs)
    lecs.sort()
    return lecs


def _pair_entropy(ecs_pair):
    """Calculates the entropy based similarity (distance) of a pair
    of EC numbers. Returns the entropy (H) of the ec pair.

    DEPREACTED: This function does NOT reflect the hierarchy nature of
    EC classification, use function pair_entropy_hierarchy() instead.

    Arguments:
    - `ecs_pair`: a pair of EC numbers (tuple), EC3 (3 levels of clasification)
    """
    f1, f2, f3 = 15., 10., 1.  # factors to ponder the value of each level
    ec1 = ecs_pair[0].split('.')
    ec2 = ecs_pair[1].split('.')
    entropies = []
    for i in range(len(ec1)):
        if ec1[i] == ec2[i]:
            entropies.append(0)  # if both elements in the column are
# equal, H = 0
        else:
            entropies.append(1)  # if both elements in the column are
# dissimilar, H = 1
    # global normalized entropy of the ecs
    H = (f1 * entropies[0] + f2 * entropies[1] +
         f3 * entropies[2]) / (f1 + f2 + f3)
    return H


def pair_entropy_hierarchy(ecs_pair):
    """Calculates the entropy based similarity (distance) of a pair
    of EC numbers.Returns the entropy (H) of the ec pair. This function 
    takes into account the hierarchy of EC classification. Therefore 
    if the first EC level do not coincide, the H=0, despite the coincidence
    of the other levels. This is the prefered function to create the 
    EC similarity matrix.

    Arguments:
    - `ecs_pair`: a pair of EC numbers (tuple), EC3 (3 levels of clasification)
    """
    f1, f2, f3 = 15., 10., 1.
    ec1 = ecs_pair[0].split('.')
    ec2 = ecs_pair[1].split('.')
    entropies = []
    if ec1[0] == ec2[0] and ec1[1] == ec2[1] and ec1[2] == ec2[2]:
        entropies = [0, 0, 0]
    elif ec1[0] == ec2[0] and ec1[1] == ec2[1]:
        entropies = [0, 0, 1]
    elif ec1[0] == ec2[0]:
        entropies = [0, 1, 1]
    else:
        entropies = [1, 1, 1]
    # global normalized entropy of the ecs
    H = (f1 * entropies[0] + f2 * entropies[1] +
         f3 * entropies[2]) / (f1 + f2 + f3)
    return H


def build_entropy_matrix(ecs_list):
    """Creates EC number entropy based similarity (distance) matrix from a 
    list of ec numbers.

    Arguments:
    - `ecs_list`: List of no redundant ec numbers
    """
    entropy = pair_entropy_hierarchy
    num_ecs = len(ecs_list)
    matrix = np.zeros((num_ecs, num_ecs))
    for i in range(num_ecs):
        for j in range(num_ecs):
            if i == j:
                matrix[i, j] = 0
                break
            ecs_pair = (ecs_list[i], ecs_list[j])
            H = entropy(ecs_pair)
            matrix[i, j] = H
            matrix[j, i] = H
    return matrix


def save_matrix(sub_matrix, ecs_list, filename, format='npz'):
    """Saves de similarity matrix into a file

    Arguments:
    - `sub_matrix`: EC entropy based substitution matrix 
    - `ecs_list`: List of EC numbers included in the matrix, 
                 index referenced un the matrix
    - `filename`: Name of file to store matrix
    - `format`: [npz] save to numpy binary file. txt: save
                to text file, tabular separated
    """
    assert format in ['npz', 'txt'], 'format= must be npz or txt'
    if format == 'npz':
        fname = filename + '.' + format
        np.savez(fname, matrix=sub_matrix, ecs=ecs_list)
    elif format == 'txt':
        header = '\t'.join(ecs_list)
        fname = filename + '.' + format
        np.savetxt(fname, sub_matrix, fmt='%.5f',
                   header=header, delimiter='\t')


def main():
    args = parser()
    ecs = load_ec_file(args.ecfile, ecl=3)
    mat = build_entropy_matrix(ecs)
    fname = 'h_ent_mat'
    save_matrix(mat, ecs, fname, format=args.store_format)
    if args.plot:
        import matplotlib.pyplot as plt
        plt.matshow(mat)
        plt.colorbar(label='Entropy')
        plt.savefig('ent_matrix.png', dpi=200)
        plt.close()


if __name__ == '__main__':
    main()
    # from getEC3fromDB import EClist
    # import sqlite3 as s3
    # from sys import argv

    # db = s3.connect(argv[1])

    # ecs = EClist(db)

    # # No hierarchy
    # mat = build_entropy_matrix(ecs, function=0)
    # # Hierarchy
    # hmat = build_entropy_matrix(ecs, function=1)

    # # sabe matrices
    # save_matrix(mat, ecs, 'entropy_matrix')
    # save_matrix(hmat, ecs, 'entropy_matrix_hierarchy')
