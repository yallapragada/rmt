'''
Created on Mar 23, 2015
@author: uday
implementation of https://www.pnas.org/content/108/28/11530.full for HA sequences in influenza
'''

import numpy as np
from sequence_tools import read_fasta
from collections import Counter
from scipy import stats
from sys import argv
import matplotlib.pyplot as plt
import pylab


def map_sites(sites, residue_dict):
    """ maps a list of site positions (after pruning) to original site positions
    """
    new_sites = []
    for site in sites:
        new_sites.append(residue_dict[site])
    return new_sites


def map_residue_positions(novars, pruned_sequence_length):
    """ adds positions of novars residues and gives original residue positions
    """
    counter = 0
    i = 0
    residue_dict = {}
    distinct_novars = list(set(novars))
    distinct_novars.sort()
    while counter < pruned_sequence_length:
        if (counter in distinct_novars):
            residue_dict[counter] = i + 1
            i = i + 2
        else:
            residue_dict[counter] = i
            i = i + 1
        counter = counter + 1
    return residue_dict


def filter_and_generate_binary(data):
    """Removes malformed strains and all sites that have non-varying residues or
    gaps.
 
    The return value x is a 2D boolean array, where:
        *each row corresponds to a residue
        *each column corresponds to a different strain
        *a 1 means that the given strain at that residue is identical to the
            consensus sequence
    """
    # Remove strains that aren't similar enough to each other
    data = filter_strains(data)

    # Find the most common residue at each nucleotide location
    consensus_sequence, novars = get_consensus(data)

    # exclude all sites that have a gap in the consensus strain
    gaps = [idx for idx, res in enumerate(consensus_sequence) if res is '-']
    novars.extend(gaps)

    data, consensus_sequence = strip_positions(data, consensus_sequence, novars)
    x = generate_binary_matrix(data, consensus_sequence)

    return x, novars


def filter_strains(data):
    """ Filters out strains that are sufficiently different from the rest
    """
    length_counter = Counter([len(strain) for strain in data])
    most_common = length_counter.most_common()[0][0]

    # Collect only strains that are the right length into the return variable
    good_data = []
    for sequence in data:
        if len(sequence) == most_common:
            good_data.append(sequence)

    return good_data


def get_consensus(strains):
    """ Get the consensus sequence of aligned strains
    """

    residue_counters = [Counter() for residue in strains[0]]

    for strain in strains:
        for index, residue in enumerate(strain):
            residue_counters[index][residue] += 1

    consensus_list = [counter.most_common()[0][0]
                      for counter in residue_counters]

    # list of positions with no variation
    novars = []
    for i, counter in enumerate(residue_counters):
        if len(counter) == 1:
            novars.append(i)

    consensus = ''.join(consensus_list)

    return consensus, novars


def strip_positions(data, consensus, novar):
    """Remove positions given in novar from all of the strains as well as the
    consensus strain.
    """
    data = [strip_single_position(seq, novar) for seq in data]
    consensus = strip_single_position(consensus, novar)

    return data, consensus


def strip_single_position(string, novar):
    "Remove positions given in novar from a single string"
    novar = set(novar)
    return "".join([char for i, char in enumerate(string)
                    if i not in novar])


def generate_binary_matrix(data, consensus):
    """ Generates a binary array x_i(s), where:
        * Each column corresponds to a particular strain
        * Each row corresponds to a particular site
        * The element is 1 if that strain at that site is indentical to the
           consensus sequence at that site
    """

    x = np.zeros((len(consensus), len(data)), dtype=bool)
    for s, strain in enumerate(data):
        for i, site in enumerate(strain):
            x[i, s] = (site == consensus[i])

    return x


def clean_matrix(correlation_matrix, lambda_cutoff):
    """ Uses RMT to clean the correlation matrix
    Every eigenvector with an eigenvalue greater than the cutoff is used to
    generate a new correlation matrix
    """
    # Calculate the eigenvalues and eigenvectors
    eigvals, vecs = np.linalg.eigh(correlation_matrix)
    plt.hist(eigvals, bins=200)  # histogram with 200 bins
    plt.ylim(0, 25)
    pylab.show()

    max_eig = np.sort(eigvals)[::-1]
    print ("top 10 eigen values of correlation matrix ", max_eig[:10])

    # create a clean matrix using RMT
    clean = np.zeros_like(correlation_matrix)
    for k, eigval in enumerate(eigvals):
        if eigval > lambda_cutoff and eigval != max(eigvals):
            # For each eigenvalue larger than the cutoff, compute the outer
            # product, and add it to the clean matrix. This is equation S5 in RMT paper
            clean += eigval * np.outer(vecs[:, k], vecs[:, k])
    return clean


def clean_phylogeny(binary_matrix):
    """ Cleans the binary matrix by removing the contribution of phylogeny
    This is section S4 of RMT paper.
    """

    eigvals, eigvecs = np.linalg.eigh(np.corrcoef(binary_matrix))

    # eigenvector corresponding to the largest eigen value
    u1 = eigvecs[:, eigvals.argmax()]

    num_residues, num_strains = np.shape(binary_matrix)

    # Equation S11
    M = np.array([sum(u1[i] * binary_matrix[i, s]
                      for i in range(num_residues))
                  for s in range(num_strains)])

    alpha = np.zeros((num_residues, num_strains))
    beta = np.zeros(num_residues)

    for i in range(num_residues):
        # "The value of the parameters alpha_i and beta_i are estimated through
        # a least square regression..."
        slope, intercept, r_value, p_value, std_err = stats.linregress(M, binary_matrix[i, :])
        alpha[i, :] = intercept
        beta[i] = slope

    # Equation S10:
    epsilon = binary_matrix - alpha - np.outer(beta, M)

    # Equation S12:
    return alpha + epsilon


def find_distinct_evo(binary_matrix):
    """ Removes evolutionarily distinct sequences
 
    Calculates a correlation matrix for the strains as they relate to each
    other, and then removes those that are significantly different
 
    This is section S5 of the paper
    """
    gamma = np.cov(binary_matrix.T)
    eigvals, vecs = np.linalg.eigh(gamma)
    vecs = vecs.T

    # Using the projections along eigenvector2 with a empirically determined cut-off
    proj1 = [np.dot(gamma[i], vecs[-1]) for i in range(len(eigvals))]
    proj2 = [np.dot(gamma[i], vecs[-2]) for i in range(len(eigvals))]
    return [pos for pos, proj in enumerate(proj2) if proj > -.1]


def find_cutoff(alignment):
    eigs = []
    alignment = alignment.copy()
    nresidues, nstrains = np.shape(alignment)

    global allcorrs
    allcorrs = np.empty(100 * nresidues ** 2)

    for i in range(100):
        # Shuffle the residues at each position
        for r in range(nresidues):
            alignment[r, :] = np.random.permutation(alignment[r, :])

        # Calculate the correlation coefficient
        corr = np.corrcoef(alignment, bias=1)

        # Add the eigenvalues to the running list of eigenvalues
        eigs.extend(np.linalg.eigvalsh(corr))

        allcorrs[i * nresidues ** 2:(i + 1) * nresidues ** 2] = \
            abs((corr * ~np.identity(nresidues, dtype=bool)).ravel())

    return eigs


def determine_sectors(correlation_matrix, lambda_cutoff):
    """ Determines the sectors of the protein
 
    See sections S6 and S7 of the supplementals.
 
    This function returns both the strongest eigenvalue at a given residue and
    the a list of counters with the projection of each residue onto significant
    eigenvectors
    """

    eigvals, vecs = np.linalg.eigh(correlation_matrix)

    n_residues = n_vectors = len(eigvals)

    loadings = [Counter() for i in range(n_residues)]

    # removes the autocorrelations, which should typically be much higher than
    # the inter-residue correlations
    # This works by multiplying by the inverse of the identity matrix
    othercorrs = abs(correlation_matrix
                     * (~ np.identity(n_residues, dtype=bool)))

    for r in range(n_residues):
        if max(othercorrs[r]) < 0.15:
            # "We chose to exclude from sectors those residues that did not
            # show any correlation higher than 0.15 (in absolute magnitude) with
            # any other sector residue"
            continue

        for k in range(n_vectors):
            if eigvals[k] > lambda_cutoff:
                loading = np.dot(correlation_matrix[:, r], vecs[:, k])
                loadings[r][k] = abs(loading)

    best = [(l.most_common()[0][0] if (len(l) > 0) else None)
            for l in loadings]

    return best, loadings


def create_sectors(best, residue_map):
    sectors = []
    for sec in set(best):
        if sec is None:
            continue

        sites = []
        for res in range(len(best)):
            if best[res] == sec:
                sites.append(res)
        sectors.append(sites)
        print (map_sites(sites, residue_map))
    return sectors


def sortSectors(loadings, sectors, corr_matrix_clean):
    bestloadings = []
    for i, l in enumerate(loadings):
        if len(l) > 0:
            bestloadings.append(l.most_common()[0][1])
        else:
            bestloadings.append(0)

    bestloadings = np.array(bestloadings)

    sorted_secs = []
    for sec in sectors:
        sec = np.array(sec)
        ls = bestloadings[sec]
        sorted_secs.append(sec[ls.argsort()])

    seclist = np.hstack(sorted_secs)

    plt.imshow(corr_matrix_clean[np.meshgrid(seclist, seclist)])
    pylab.show()


if __name__ == '__main__':
    flu_seq_file = argv[1]
    flu_fasta = read_fasta(flu_seq_file)
    flu_data = [flu_fasta[name] for name in flu_fasta]

    # convert to binary
    strains, novars1 = filter_and_generate_binary(flu_data)

    # remove duplicates
    distinct_strains = find_distinct_evo(strains)
    flu_data2 = [strain for idx, strain in enumerate(flu_data) if idx in
                 distinct_strains]

    # filter novars
    x, novars2 = filter_and_generate_binary(flu_data2)

    # create a residue map
    residue_map = map_residue_positions(novars2, rows)

    # remove phylogeny
    x = clean_phylogeny(x)

    # correlation matrix and eigen values
    corr_matrix = np.corrcoef(x, bias=1)
    eigs = find_cutoff(x)
    lambda_cutoff = max(eigs)

    # clean matrix using RMT
    corr_matrix_clean = clean_matrix(corr_matrix, 4.0)

    # determine sectors
    best, loadings = determine_sectors(corr_matrix_clean, lambda_cutoff)
    sectors = create_sectors(best, residue_map)
    sortSectors(loadings, sectors, corr_matrix_clean)
