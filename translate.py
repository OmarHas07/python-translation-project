#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the translated amino acids.
    """
    rna_sequence = rna_sequence.upper()
    amino_acids = []
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        if codon in genetic_code:
            amino_acid = genetic_code[codon]
            if amino_acid == '*':
                break
            amino_acids.append(amino_acid)
    return ''.join(amino_acids)
"""
the above code I got from google, it was a little bit complicated for me but I'm now understand it, because the nucleotides are 
lowercase rna_seq. upper is used to make sure the all the input is upper case, creating an empty variable for the AAs, 
then in the for loop, (0 to start from the very begining, across the length of rnaseq, and in 3 steps of bases. 
that defines the codon to start from first AA, and take 3 Bases or 3 steps to take 3 bases for each codon then to get the AA from the codon in the gentic code
which already defined previously, if AA is a stop codon, then it breaks, 

"""

def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    list
        A list of strings; each string is an sequence of amino acids encoded by
        `rna_sequence`.
    """

    rna_sequence = rna_sequence.upper()
    all_translations = []

    for frame in range(3):
        amino_acids = []
        started_translation = False
        for i in range(frame, len(rna_sequence), 3):
            codon = rna_sequence[i:i+3]
            if codon in genetic_code:
                amino_acid = genetic_code[codon]
                
                if not started_translation and amino_acid == 'M':
                    started_translation = True
                    
                if started_translation:
                    if amino_acid == '*':
                        break
                    amino_acids.append(amino_acid)
        
        if started_translation:
            all_translations.append(''.join(amino_acids))

    return all_translations




def get_reverse(sequence):
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_reverse('AUGC')
    'CGUA'
    """
    rev = sequence[::-1]
    return(rev).upper()


def get_complement(sequence):
    """Get the complement of a `sequence` of nucleotides.

    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_complement('AUGC')
    'UACG'
    """
    sequence = sequence.upper() 
    complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    complement = [complement_map[base] for base in sequence]
    return ''.join(complement)
"""
the first sequence just ensure that the input is upper case, comp map, dicatianary to the RNA genetic code, then using the 
dict I created (comp_map) to go through the sequence and get the complement, last step is to convert the list into a single 
string with no separator between elements. 
"""

def reverse_and_complement(sequence):
    """Get the reversed and complemented form of a `sequence` of nucleotides.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> reverse_and_complement('AUGC')
    'GCAU'
    """
    sequence = sequence.upper() 
    complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    complement = [complement_map[base] for base in sequence]
    reverse_complement = ''.join(complement[::-1])
    return reverse_complement

def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the longest sequence of amino acids encoded by
        `rna_sequence`.
    """





"""
test the six reading frames which is the forward the reverse geting the longest sequence maybe using leng function     

"""

def get_longest_peptide(rna_sequence, genetic_code):
    def get_peptides(sequence):
        peptides = []
        for frame in range(3):
            peptide = []
            started_translation = False
            for i in range(frame, len(sequence), 3):
                codon = sequence[i:i+3]
                if codon in genetic_code:
                    amino_acid = genetic_code[codon]
                    if not started_translation and amino_acid == 'M':
                        started_translation = True
                    if started_translation:
                        if amino_acid == '*':
                            peptides.append(''.join(peptide))
                            peptide = []
                            started_translation = False
                        else:
                            peptide.append(amino_acid)
                            
            # Add the last peptide if the sequence did not end with a stop codon.
            if started_translation and peptide:
                peptides.append(''.join(peptide))
        return peptides

##I couldn't use the reverse complenment function I already made earlier it kept giving me an error for some reason, I thought I can just use it instead of writing a new function here
    def reverse_complement(sequence):
        complement = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement[base] for base in reversed(sequence))

    rna_sequence = rna_sequence.upper()
    reverse_complement_rna = reverse_complement(rna_sequence)

    # Get peptides for both forward and reverse reading frames.
    peptides = get_peptides(rna_sequence) + get_peptides(reverse_complement_rna)

    # Find the longest peptide among all peptides.
    longest_peptide = max(peptides, key=len, default='')

    return longest_peptide




if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
