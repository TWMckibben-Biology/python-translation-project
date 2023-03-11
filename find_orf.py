#! /usr/bin/env python3

import sys
import re

def vet_nucleotide_sequence(sequence):
    #REGEX specific to either RNA or DNA
    rna_pattern_str = r'[AUCGaucg]*$'
    dna_pattern_str = r'[ATCGatcg]*$'
    
    #Sets the REGEX str as a pattern
    rna_pattern = re.compile(rna_pattern_str)
    dna_pattern = re.compile(dna_pattern_str)

    #Matches the called sequence to the pattern for RNA and DNA, calls an exception if it doesnt match
    if rna_pattern.match(sequence):
        return
    if dna_pattern.match(sequence):
        return
    else:
        raise Exception("Invalid sequence: {0!r}".format(sequence))



def vet_codon(codon):
    
    #REGEX that calls three bases of RNA
    codon_pattern_str = r'[AUGCaugc][AUGCaugc][AUGCaugc]$'
    
    #Sets REGEX str as pattern
    codon_pattern = re.compile(codon_pattern_str)

    #Checks that the called sequence has codons matching RNA, returns an exception otherwise
    if codon_pattern.match(codon):
        return
    else:
        raise Exception("Invalid codon: {0!r}".format(codon))


def find_first_orf(sequence,
        start_codons = ['AUG'],
        stop_codons = ['UAA', 'UAG', 'UGA']):
    
    ###Start of our first real script

    #Make sure the sequence is valid
    vet_nucleotide_sequence(sequence)

    #Make sure the codons are valid
    for codon in start_codons:
        vet_codon(codon)
    for codon in stop_codons:
        vet_codon(codon)

    #Sets sequences to upper case
    seq = sequence.upper()
    starts = [c.upper() for c in start_codons]
    stops = [c.upper() for c in stop_codons]
    
    #Double checks that our DNA has been converted to RNA
    seq = seq.replace('T', 'U')

    #Checks for a start codon, looks at the codons after that, finishes at a stop codon
    orf_pattern_str = r'(' + '|'.join(start_codons) + ')([GUAC]{3})*(' + '|'.join(stop_codons) +')'
    
    #Sets REGEX str as pattern
    orf_pattern = re.compile(orf_pattern_str)
    
    #Searches the sequence using orf_pattern and saves it as a variable, returns it in a string
    match_object = orf_pattern.search(seq)
    if match_object:
        return match_object.group()
    return ''


def parse_sequence_from_path(path):
    
    #Checks that file can be found by given path, throws an exception otherwise
    try:
        file_stream = open(path, 'r')
    except FileNotFoundError as e:
        sys.stderr.write("Sorry, couldn't find path {}".format(path))
        raise e
    except IsADirectoryError as e:
        sys.stderr.write("Sorry, path {} appears to be a directory".format(
                path))
        raise e
    except:
        sys.stderr.write("Sorry, something went wrong when trying to open {}".format(
                path))
        raise
    #Saves the match_object.group output as the variable "sequence"
    sequence = ''
    
    #Basically removes whitespace that would break downstream code
    for line in file_stream:
        sequence += line.strip()
    return sequence


def main():
    import argparse

    ###Here is the meat and potatoes, this is the compilation of all previous scripts
    # Create a command-line parser object
    parser = argparse.ArgumentParser()

    #Set the start and stop codons
    default_start_codons = ['AUG']
    default_stop_codons = ['UAA', 'UAG', 'UGA']

    # Tell the parser what command-line arguments this script can receive
    parser.add_argument('sequence',
            metavar = 'SEQUENCE',
            type = str,
            help = ('The sequence to search for an open-reading frame. '
                    'If the path flag (\'-p\'/\'--path\') is specified, '
                    'then this should be a path to a file containing the '
                    'sequence to be searched.'))
    parser.add_argument('-p', '--path',
            action = 'store_true',
            help = ('The sequence argument should be treated as a path to a '
                    'containing the sequence to be searched.'))
    parser.add_argument('-s', '--start-codon',
            type = str,
            action = 'append', # append each argument to a list
            default = None,
            help = ('A start codon. This option can be used multiple times '
                    'if there are multiple start codons. '
                    'Default: {0}.'.format(" ".join(default_start_codons))))
    parser.add_argument('-x', '--stop-codon',
            type = str,
            action = 'append', # append each argument to a list
            default = None,
            help = ('A stop codon. This option can be used multiple times '
                    'if there are multiple stop codons. '
                    'Default: {0}.'.format(" ".join(default_stop_codons))))

    # Parse the command-line arguments into a 'dict'-like container
    args = parser.parse_args()

    # Check to see if the path option was set to True by the caller. If so, parse
    # the sequence from the path
    if args.path:
        sequence = parse_sequence_from_path(args.sequence)
    else:
        sequence = args.sequence

    # Check to see if start/stop codons were provided by the caller. If not,
    # use the defaults.
    if not args.start_codon:
        args.start_codon = default_start_codons
    if not args.stop_codon:
        args.stop_codon = default_stop_codons

    orf = find_first_orf(sequence = sequence,
            start_codons = args.start_codon,
            stop_codons = args.stop_codon)
    sys.stdout.write('{}\n'.format(orf))


if __name__ == '__main__':
    main()
