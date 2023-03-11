#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    
    #Counts the length of the sequence, seperate it into codons, reference it against the dictionary, return the protein translation
    rna_sequence = rna_sequence.upper()
    protein = "" 
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        codon = codon.upper()
        if len(codon) < 3 :
            break
        if genetic_code[codon] == "*" :
            break
        protein += genetic_code[codon]   
    return protein 
    
def get_all_translations(rna_sequence, genetic_code):
    
    #Sets a loop that progresses through each base and gets all possible translations
    rna_sequence = rna_sequence.upper()
    print(rna_sequence)
    aa = ""
    protein_list = []
    protein = ""
    start = ""
    start_list = []

    for i in range(0, len(rna_sequence), 1):
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
    print("protein is...", protein)            
    print(protein_list)
    return(protein_list)
    
def get_reverse(sequence):
    
    #Returns an uppercase, reversed sequence
    stringlength=len(sequence)
    sequence=sequence[stringlength::-1]
    sequence=sequence.upper()
    return(sequence)

def get_complement(sequence):
    
    #Replaces A's with U's, etc, makes the complement of 'sequence'
    sequence=sequence.upper()
    complement = {"A":"U", "U":"A", "C":"G", "G":"C"}
    bases = list(sequence)
    bases = [complement[base] for base in bases]
    print("".join(bases))
    return "".join(bases)

def reverse_and_complement(sequence):
    
    #Combines previous code to reverse 'sequence' and then get its complement
    stringlength=len(sequence)
    sequence=sequence[stringlength::-1]
    sequence=sequence.upper()
    complement = {"A":"U", "U":"A", "C":"G", "G":"C"}
    bases = list(sequence)
    bases = [complement[base] for base in bases]
    print("".join(bases))
    return "".join(bases)

def get_longest_peptide(rna_sequence, genetic_code):
    
    #Checks all rna sequences to see which ones make the longest valid proteins, must start with a start codon
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

    for i in range(0, len(rna_sequence), 1):
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
    
    reverse=reverse_and_complement(rna_sequence)
    
    for a in range(0, len(reverse), 1):
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
