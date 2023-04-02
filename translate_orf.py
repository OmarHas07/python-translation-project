
import sys 
import re  
#import vet_codon
#import find_first_orf 
#import parse_sequence_from_path

def vet_nucleotide_sequence(sequence):


    rna_pattern_str = r'[AUCGaucg]*$'
    dna_pattern_str = r'[ATCGatcg]*$'

    rna_pattern = re.compile(rna_pattern_str)
    dna_pattern = re.compile(dna_pattern_str)

    if not sequence.strip(): # empty or whitespace-only string is valid
        return

    if not rna_pattern.match(sequence) and not dna_pattern.match(sequence):
        raise Exception("Invalid sequence: {0!r}".format(sequence))

    if rna_pattern.match(sequence) and dna_pattern.match(sequence):
        raise Exception("Invalid sequence: {0!r}: sequence contains both RNA and DNA bases".format(sequence))

    if rna_pattern.match(sequence):
        return
    elif dna_pattern.match(sequence):
        return

def vet_codon(codon):
    codon_pattern_str = r'^[AUGC]{3}$'

    codon_pattern = re.compile(codon_pattern_str, re.IGNORECASE)

    if codon_pattern.match(codon):
        return
    else:
        raise Exception("Invalid codon: {0!r}".format(codon))


def find_first_orf(sequence,
        start_codons = ['AUG'],
        stop_codons = ['UAA', 'UAG', 'UGA']):
    
    # Make sure the sequence is valid
    vet_nucleotide_sequence(sequence)

    # Make sure the codons are valid
    for codon in start_codons:
        vet_codon(codon)
    for codon in stop_codons:
        vet_codon(codon)

    # Get copies of everything in uppercase
    seq = sequence.upper()
    starts = [c.upper() for c in start_codons]
    stops = [c.upper() for c in stop_codons]
    # Make sure seq is RNA
    seq = seq.replace('T', 'U')

    orf_pattern_str = r'AUGGUAUAA'


    orf_pattern_str = r'(' + '|'.join(starts) + r')([ACGU]{3})*(' + '|'.join(stops) + r')'
    orf_pattern = re.compile(orf_pattern_str)

    # Search the sequence
    match_object = orf_pattern.search(seq)
    if match_object:
        return match_object.group()
    return ''





def translate_sequence(sequence):
    # Using dict to define the aa each codon I don't know if there is a smarter way of doing it. 
    codon_table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    # Make sure the sequence is valid
    vet_nucleotide_sequence(sequence)

    # Make sure the sequence has a valid ORF
    orf = find_first_orf(sequence)
    if not orf:
        raise Exception("Invalid sequence: {0!r}: no valid ORF found".format(sequence))

    # Translate the sequence
    aa_sequence = ''
    for i in range(0, len(orf), 3):
        codon = orf[i:i+3]
        aa = codon_table.get(codon, 'X')
        aa_sequence += aa

    return aa_sequence

