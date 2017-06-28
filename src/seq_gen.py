#!/usr/bin/env python2
# Generator for binary sequences using MSAs
# Written by Arjun Srinivasan

from numpy import array,transpose,savetxt
from scipy.stats import mode
import sys

from Bio import AlignIO

def generate_binary_sequence(fasta_list):
    """ Generates a binary sequence out of an MSA.

    Arguments:
    fasta_list -- A list of protein sequences.
    """
    mod = [mode(x)[0][0] for x in transpose(array([list(z) for z in fasta_list]))]
    return array([[1 if x==mod[i] else 0 for i,x in enumerate(y)] for y in fasta_list])

# test code below here
if __name__ == '__main__':
    bin_seq = [x.seq for x in AlignIO.read(sys.argv[1], 'fasta')]
    savetxt(sys.argv[1][:-6]+"_bin.out", bin_seq)
    print(generate_binary_sequence(bin_seq))