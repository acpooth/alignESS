#!/usr/bin/python
# 
#
# @uthor: acph
#  


def timing():
    scores = []
    aligns = []
    t1 = time()
    for i in range(len(seqs)):
        s1 = seqs[i]
        for j in range(i):
            s2 = seqs[j]
            result = NW.NW(mat, ecs, s1,s2)
            scores.append(result[2])
            aligns.append((result[0], result[1]))
    elapsed = time() - t1
    return scores, aligns, elapsed

def timingh():
    scores = []
    aligns = []
    t1 = time()
    for i in range(len(seqs)):
        s1 = seqs[i]
        for j in range(i):
            s2 = seqs[j]
            result = NW.NW(hmat, ecs, s1,s2)
            scores.append(result[2])
            aligns.append((result[0], result[1]))
    elapsed = time() - t1
    return scores, aligns, elapsed

