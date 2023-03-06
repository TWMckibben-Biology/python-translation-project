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
    #print(rna_sequence)
   # rnalen = length(rna_sequence)
    rna_sequence = rna_sequence.upper()
    protein = "" 
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        #print(codon)
        codon = codon.upper()
        if len(codon) < 3 :
            break
        if genetic_code[codon] == "*" :
            #print("identified a stop codon")
            break
        protein += genetic_code[codon]   
    
    #print(protein)
            #protein = protein.upper()
    return protein 

    #print(rna_bases)
    #rna_translated = "".join(rna_bases)
    #print(rna_translated)
    #return "".join(rna_bases)
    
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
    print(rna_sequence)
    aa = ""
    protein_list = []
    protein = ""
    start = ""
    start_list = []
#    import pdb; pdb.set_trace()
    for i in range(0, len(rna_sequence), 1):
        #frame1 = rna_sequence[i:i+3]
        if rna_sequence[i:i+3] == "AUG" :
            start_list.append(i)
    for y in start_list:
        for x in range(y, len(rna_sequence), 3):
            codon = rna_sequence[x:x+3]
            print(codon + "-")
            if len(codon) < 3:
                print("codon less than 3 identified")
                break        
            if genetic_code[codon] == "*" :
                print("identified a stop codon")
                break
            protein += genetic_code[codon]
        protein_list.append(protein)
        protein = ""
#    for z in range(len(rna_sequence) - 2):
#        codon = rna_sequence[z:z+3]
#        if codon == "AUG" :
#            remaining_seq = rna_sequence[z:len(rna_sequence)]
#            for x in range(0, len(remaining_seq), 3):
#                aa = genetic_code[remaining_seq]
#            if aa:
#                protein_list.append(aa)
#                print(aa)
    print("protein is...", protein)            
        #break
    
    print(protein_list)
    return(protein_list)
    
def get_reverse(sequence):
    """Reverse orientation of `sequence`.
    Returns a string with `sequence` in the reverse order.
    If `sequence` is empty, an empty string is returned.
    Examples
    --------
    >>> get_reverse('AUGC')
n_list = []
    'CGUA'
    """
    #print("\n"+sequence)
    #type(sequence)
    stringlength=len(sequence)
    sequence=sequence[stringlength::-1]
    sequence=sequence.upper()
    return(sequence)

def get_complement(sequence):
    """Get the complement of a `sequence` of nucleotides.
    Returns a string with the complementary sequence of `sequence`.
    If `sequence` is empty, an empty string is returned.
    Examples
    --------
    >>> get_complement('AUGC')
    'UACG'
    """
    sequence=sequence.upper()
    complement = {"A":"U", "U":"A", "C":"G", "G":"C"}
    #print(complement[sequence])
    bases = list(sequence)
    bases = [complement[base] for base in bases]
    print("".join(bases))
    return "".join(bases)

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
    stringlength=len(sequence)
    sequence=sequence[stringlength::-1]
    sequence=sequence.upper()
    complement = {"A":"U", "U":"A", "C":"G", "G":"C"}
    bases = list(sequence)
    bases = [complement[base] for base in bases]
    print("".join(bases))
    return "".join(bases)

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
    rna_sequence = rna_sequence.upper()
    print(rna_sequence)
    longestpeptide = ""
    aa = ""
    protein_list = []
    protein = ""
    revprotein = ""
    start = ""
    revstart = ""
    rev_start_list = []
    start_list = []
#    import pdb; pdb.set_trace()
    for i in range(0, len(rna_sequence), 1):
        #frame1 = rna_sequence[i:i+3]
        if rna_sequence[i:i+3] == "AUG" :
            start_list.append(i)
    for y in start_list:
        for x in range(y, len(rna_sequence), 3):
            codon = rna_sequence[x:x+3]
            print(codon + "-")
            if len(codon) < 3:
                print("codon less than 3 identified")
                break
            if genetic_code[codon] == "*" :
                print("identified a stop codon")
                break
            protein += genetic_code[codon]
        protein_list.append(protein)
        protein = ""
    
    #stringlength=len(rna_sequence)
    #reverse=rna_sequence[stringlength::-1]
    #print(rna_sequence)
    #print(reverse)
    reverse=reverse_and_complement(rna_sequence)
    
    for a in range(0, len(reverse), 1):
        #frame1 = rna_sequence[i:i+3]
        if reverse[a:a+3] == "AUG" :
            rev_start_list.append(a)
    for b in rev_start_list:
        for c in range(b, len(reverse), 3):
            revcodon = reverse[c:c+3]
            print(revcodon + "-")
            if len(revcodon) < 3:
                print("codon less than 3 identified")
                break
            if genetic_code[revcodon] == "*" :
                print("identified a stop codon")
                break
            revprotein += genetic_code[revcodon]
        protein_list.append(revprotein)
        revprotein = ""


    print("protein list is...", protein_list)
    if protein_list:
        longestpeptide = max(protein_list, key=len)
    print(longestpeptide)
    return longestpeptide


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
