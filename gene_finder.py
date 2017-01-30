# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        return 'not a nucleotide'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverse_compelment = ''
    for i in range(1, len(dna)+1):
        # iterates backwards through the loop
        new_i = len(dna) - i
        # gets nucleotide at end of list
        nucleotide = dna[new_i]
        # finds complement of nucleotide
        complement = get_complement(nucleotide)
        # adds complement to return value reverse_compelment
        reverse_compelment = reverse_compelment + complement
    return reverse_compelment


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # Sets default ORF to start codon
    ORF = dna[0:3]
    # Defines length codons as number of full codons
    length_codons = len(dna)//3
    # iterates through codons starting with second
    for i in range(2, length_codons+1):
        # defines stant and end of condon in dna
        end = 3*i
        start = end - 3
        # gets codon from dna
        codon = dna[start:end]
        # checks for stop codon
        if codon == 'TGA' or codon == 'TAG' or codon == 'TAA':
            return ORF
        # adds codon to ORF
        else:
            ORF = ORF + codon
            # if current codon is the last, add the rest of dna strand
            # to ORF
            if (i*3 + 3) > len(dna):
                ORF = ORF + dna[end:len(dna)]
    return ORF


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    # initilizes running_length to keep track of index in dna
    running_length = 0
    # initilizes return variable list_ORFs
    list_ORFs = []
    # loops through dna sequence
    while running_length <= len(dna):
        # start where last ORF ended
        start = running_length
        # end of codon is three away from start
        end = start + 3
        # pulls codon from dna
        codon = dna[start:end]
        # checks for start codon
        if codon == 'ATG':
            # get new ORF
            new_ORF = rest_of_ORF(dna[running_length:len(dna)])
            # update return variable list_ORFs
            list_ORFs = list_ORFs + [new_ORF]
            # updates running_length to index after last ORF
            running_length = running_length + len(new_ORF)
        else:
            # updates running_length to index of next codon
            running_length = running_length + 3
    return list_ORFs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    frame_1 = dna
    frame_1_list = find_all_ORFs_oneframe(frame_1)
    frame_2 = dna[1:len(dna)]
    frame_2_list = find_all_ORFs_oneframe(frame_2)
    frame_3 = dna[2:len(dna)]
    frame_3_list = find_all_ORFs_oneframe(frame_3)
    list_ORFs = frame_1_list + frame_2_list + frame_3_list
    return list_ORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    r_compelment = get_reverse_complement(dna)
    dna_ORFs = find_all_ORFs(dna)
    reverse_compelment_ORFs = find_all_ORFs(r_compelment)
    all_ORFs = dna_ORFs + reverse_compelment_ORFs
    return all_ORFs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    rest_of_ORF("ATGAGATAGG")
