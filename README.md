# alignESS.

A program for Enzymatic Step Sequence (ESS) alignment using the Dynamic Programming (DP) Needleman–Wunsch algorithm.

The program can perform the following alignments:

    pair                ESS command line pairwise comparisson using DP
    dbalign             ESSs database(s) alignment using DP
    multi               ESSs multiple alignment using GA

The pairwise alignments are generated with the DP algorithm, and the multiple ESS alignement is generated using a Genetic Algorithm (GA).

## Dependencies

- Pair-wise alignment (simple pairwise or database alignment):
 - cython
 - numpy
 - matplotlib (only for ecs_entropy.py script)

- Multiple alignment:
 - C boost library: only for compiling the alignment algorithm*
 
* This repo contains a compiled (linux 64 bit) copy of the multiple alingment algorithm that may work fine in
 the majority of linux systems. The source code of this part of the program is not yet included.

## Usage.

### Input.


__UNDER CONSTRUCTION__


## Papers.
The programs presented here were used all or in parts in the following papers.
