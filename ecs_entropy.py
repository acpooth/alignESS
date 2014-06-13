#!/usr/bin/python
# 
#
# @uthor: acph
#  

from numpy import zeros

"""
"""

def pair_entropy(ecs_pair):
    """Obtains the entropy of a tuple of ecs. 
    
    Returns the entropy (H) of the ec pair

    Arguments:
    - `ecs_pair`: a pair of ec numbers (tuple), EC3 (3 levels of clasification)
    """
    f1, f2, f3 = 15., 10., 1. # factors to ponder the value of each level
    ec1 = ecs_pair[0].split('.')
    ec2 = ecs_pair[1].split('.')
    entropies = []
    for i in range(len(ec1)):
        if ec1[i] == ec2[i]:
            entropies.append(0) # if both elements in the column are
# equal, H = 0
        else:
            entropies.append(1) # if both elements in the column are
# dissimilar, H = 1
    H = (f1*entropies[0] + f2*entropies[1] + f3*entropies[2])/(f1+f2+f3) # global normalized entropy of the ecs
    return H



def pair_entropy_hierarchy(ecs_pair):
    """Obtains the entropy of a tuple of ecs. This function takes into acount the hierarchy of ec numbers. Therefore if the first EC level don't coincide, the entropy = 0 despite the coincidence of the other levels
    
    Returns the entropy (H) of the ec pair.

    Arguments:
    - `ecs_pair`: a pair of ec numbers (tuple), EC3 (3 levels of clasification)
    """
    f1, f2, f3 = 15., 10., 1.
    ec1 = ecs_pair[0].split('.')
    ec2 = ecs_pair[1].split('.')
    entropies = []
    if ec1[0] == ec2[0] and ec1[1] == ec2[1] and ec1[2] == ec2[2]:
        entropies = [0,0,0]
    elif ec1[0] == ec2[0] and ec1[1] == ec2[1]:
        entropies = [0,0,1]
    elif ec1[0] == ec2[0]:
        entropies = [0,1,1]
    else:
        entropies = [1,1,1]
    H = (f1*entropies[0] + f2*entropies[1] + f3*entropies[2])/(f1+f2+f3) # global normalized entropy of the ecs
    return H


def build_entropy_matrix(ecs_list, function=0):
    """Creates EC a substitution matrix from a list of ec numbers, using one of to posible functions (pair_entropy or pair_entropy_hierachy)
    
    Arguments:
    - `ecs_list`: List of no redundant ec numbers
    - `function`: function to use. 0 = pair_entropy ; 1 = pair_entropy_hierarchy
    """
    if function == 0:
        entropy = pair_entropy
    elif function == 1:
        entropy = pair_entropy_hierarchy
    
    num_ecs = len(ecs_list)
    matrix = zeros((num_ecs, num_ecs))
    for i in range(num_ecs):
        for j in range(num_ecs):
            if i == j:
                matrix[i,j] = 0
                break
            ecs_pair = (ecs_list[i], ecs_list[j])
            H = entropy(ecs_pair)
            matrix[i,j] = H
            matrix[j,i] = H
            
    return matrix
    
def save_matrix(sub_matrix, ecs_list, filename):
    """Saves de substitution matrix into a numpy binary file (.npz)
    
    Arguments:
    - `sub_matrix`: EC entropy based substitution matrix 
    - `ecs_list`: List of EC numbers included in the matrix, index referenced un the matrix
    - `filename`: Name of file to store matrix
    """
    from numpy import savez
    savez(filename, matrix=sub_matrix, ecs=ecs_list)


if __name__ == '__main__':
    """  Creates the files containing the Entropy-based ECs substitution matrices. One that takes into account the hierarchy of EC numbers and other that do not take this considerations

USE: $python ecs_entropy.py [sqlite3 database file]
    
"""
    from getEC3fromDB import EClist
    import sqlite3 as s3
    from sys import argv
    

    db = s3.connect(argv[1])

    ecs = EClist(db)

    # No hierarchy
    mat = build_entropy_matrix(ecs, function=0)
    # Hierarchy
    hmat = build_entropy_matrix(ecs, function=1)

    #sabe matrices
    save_matrix(mat, ecs, 'entropy_matrix')
    save_matrix(hmat, ecs, 'entropy_matrix_hierarchy')
