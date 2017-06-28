#!/usr/bin/env python2
# Random Matrix Theory (for protein sequence analysis)
# Written by Arjun Srinivasan

from numpy import argsort,array,arange,amax,corrcoef,inner,outer,real,transpose,sum as nsum,sort as nsort,sqrt,savetxt,mean,var
from numpy.random import permutation,randint
from numpy.linalg import eig,eigvals

from seq_gen import generate_binary_sequence
from Bio import AlignIO
import sys
import multiprocessing as mp


def random_matrix_generation(tseq, num_permutations):
    """ Generates eigenvalues for a random matrix using permutations.

    Arguments:
    tseq -- Transposed list of binary sequences. Make sure this is transposed.
    num_permutations -- Number of permutations.

    """
    return array([eigvals(corrcoef(array([permutation(x) for x in tseq]))) for y in arange(num_permutations)]).flatten()

def generate_eigenspectrum(tseq):
    """ Generates the eigenspectrum.
    
    Arguments:
    tseq -- Transposed list of binary sequences.
    """
    return eig(corrcoef(tseq))

def remove_noise(evals, evecs, max_eig, phylogeny=True):
    """ Removes noise, and if necessary, phylogeny from the set of eigenvalues/eigenvectors.

    Arguments:
    evals -- Eigenvalues of correlation matrix.
    evecs -- Eigenvectors of correlation matrix.
    phylogeny -- True if phylogeny is to be removed.
    """
    return zip(*[(x,y) for x,y in zip(evals, evecs) if x > max_eig and (x != amax(evals) if phylogeny else True)])


def clean_correlation_matrix(evals, evecs, max_eig, phylogeny=True):
    """ Cleans the correlation matrix of noise, and provides an option to remove the largest eigenvector to clean the matrix of phylogeny.

    Arguments:
    evals -- Eigenvalues of correlation matrix.
    evecs -- Eigenvectors of correlation matrix.
    max_eig -- Theoretical random maximum eigenvalue. We ignore anything below the minimum eigenvalue.
    """
    return real(nsum((x*outer(y,y) for x,y in zip(evals, evecs))))



def generate_sectors(tseq, eigvectors):
    """ Generate projections of sequence vectors onto important eigenvectors.

    Arguments:
    tseq -- Transposed list of sequences.
    eigvectors -- List of important eigenvectors.
    """
    return array([[inner(x,y) for x in eigvectors] for y in tseq])

def theoretical_eig(l, n):
    """ Returns the minimum and maximum eigenvalues in the eigenspectrum for a matrix given length l and n vectors.
    """
    q = float(l)/n 
    return 1+1/q-2*sqrt(1/q), 1+1/q+2*sqrt(1/q)

def remove_zero_variance(tseq):
    """ Removes sequences with zero variance.

    Arguments:
    tseq -- Transposed list of sequences.
    """
    return array([x for x in tseq if var(x) != 0])


# test code below here
if __name__ == '__main__':
    a = [str(x.seq) for x in AlignIO.read(sys.argv[1], "fasta")]
    min_eig, max_eig = theoretical_eig(len(a),len(a[0]))
    print(max_eig)
    matrix = generate_binary_sequence(a)
    #matrix = array([[randint(0,2) for x in arange(500)] for z in arange(1600)])
    transm = remove_zero_variance(transpose(matrix))
    max_eig = nsort(random_matrix_generation(transm,2))[-1]
    print(max_eig)
    evals, evecs = generate_eigenspectrum(transm)
    cor = corrcoef(transm)
    print(len(transpose(evecs)[0]))
#    eval_imp, evec_imp = remove_noise(evals, transpose(evecs), max_eig)
    indices = argsort(evals)
    indices = indices[::-1]
    evecs = transpose(evecs[:,indices])
    evals = evals[indices]
    print(len(evecs))
#    print(len(eval_imp))
#    print(len(evec_imp[0]))
#    print(len(transm[0]))
    savetxt(sys.argv[1][:-6]+"_eigs.out", evals)
    savetxt(sys.argv[1][:-6]+"_eigvecs.out", evecs)
    savetxt(sys.argv[1][:-6]+"_projections.out", generate_sectors(cor, evecs[:6]))
