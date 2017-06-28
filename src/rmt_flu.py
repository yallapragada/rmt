#!/usr/bin/env python2
# Random Matrix Theory (for protein sequence analysis)
# singular value decomposition

from numpy import zeros,dot,delete,argsort,array,arange,amax,corrcoef,inner,outer,real,transpose,sum as nsum,sort as nsort,sqrt,savetxt,mean,var
from numpy.random import permutation,randint
from numpy.linalg import eig
from numpy.linalg import eigvals

from seq_gen import generate_binary_sequence
from Bio import AlignIO
import sys

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


def createSectors(evals, evecs, cor, max_eig):
    
    indices = argsort(evals)
    for i in range(evals.shape[0]):
        if (evals[i]>max_eig):
            print ("significant eigenval ", i, " ", indices[i], " ", evals[i])
    cor = delete(cor, 353, 0)
    cor = delete(cor, 353, 1)
    bad_vals = range(12,373)
    evecs = delete(evecs,0,1)
    evecs = delete(evecs,0,0)
    evecs = delete(evecs,bad_vals,1)
    print ("dim of new eigenvec after removing phylogeny ", evecs.shape)
    reduced = zeros( (373,8), float)
    for i in range(373):
        for j in range(8):
            reduced[i][j] = dot(cor[i], evecs[:,j].real)

    reducedT = transpose(reduced)
    
    for i in range(8):
        print (i)
        for j in range(373):
                if (reducedT[i][j] > 0.2):
                    print (j)
        print (" ")    

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
    print ("num sequences, sequence length ", len(a), len(a[0]))

    min_eig, max_eig = theoretical_eig(len(a),len(a[0]))
    print ("theoretical min and max eigen values ", min_eig, max_eig)

    matrix = generate_binary_sequence(a)
    transm = remove_zero_variance(transpose(matrix))
    print ("dimensions of input data after removing zero covariance", transm.shape)

    cor = corrcoef(transm)
    print ("dimensions of correlation matrix ", cor.shape)
    
    cor2 = corrcoef(transpose(matrix))
    print ("dimensions of initial correlation matrix ", cor2.shape)

    rand_eig_vals = random_matrix_generation(transm,50)
    print ("total no. of random eig values ", len(rand_eig_vals))
    
    max_eig = nsort(rand_eig_vals)[::-1]
    print ("top 10 eigen values of random matrix ", max_eig[:10])

    evals, evecs = generate_eigenspectrum(transm)
    indexes = argsort(evals)[::-1]
    print (indexes)
    eigval = evals[indexes]
    print (eigval[:10])
    #createSectors(evals, evecs, cor, max_eig)