#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Identify genes in a procaryote assembly or contig.
   Simple approch based on open reading frame detection and Shine-
   Dalgarno motif.
"""

import argparse
import sys
import os
import csv
import re


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='genome_file', type=isfile, required=True,
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', dest='min_gene_len', type=int,
                        default=50, help="Minimum gene length to consider")
    parser.add_argument('-s', dest='max_shine_dalgarno_distance', type=int,
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif")
    parser.add_argument('-d', dest='min_gap', type=int, default=40,
                        help="Minimum gap between two genes (shine box not included).")
    parser.add_argument('-p', dest='predicted_genes_file', type=str,
                        default=os.curdir + os.sep +"predict_genes.csv",
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=str,
                        default=os.curdir + os.sep + "genes.fna",
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file):
    """Extract the complete genome sequence as a single string.
    If several only one sequence is considered
    """
    sequence = ""
    getone = False
    with open(fasta_file, 'rt') as my_file:
        for line in my_file:
            if line.startswith(">"):
                getone = True
            else:
                sequence += line.strip().upper()
    if getone:
        return sequence
    sys.exit("No sequence found")


def find_start(start_regex, sequence, start, stop):
    """Find the start codon.
    Regex.search() return None if no match found.
    """
    found = start_regex.search(sequence, start, stop)
    if found is None:
        return found
    return found.start(0)


def find_stop(stop_regex, sequence, start):
    """Find the stop codon in same Open Reading Frame.
    """
    start_frame = start % 3
    matches = stop_regex.finditer(sequence, start)
    for match in matches:
        position = match.start(0)
        match_frame = position % 3
        if match_frame == start_frame:
            return position
    return None


def has_shine_dalgarno(shine_regex, sequence,
                       start, max_shine_dalgarno_distance, verbose = False):
    """Find a shine dalgarno motif before the start codon.
    """
    window_start = max(start - max_shine_dalgarno_distance, 0)
    if verbose:
        print(f"search will begin at position {window_start}")
    options = shine_regex.finditer(sequence, window_start, start)
    for option in options:
        if verbose:
            print(f"Current tested position is {option}")
            print(f"match ending position: {option.end()}")
            print(f"start - end option position: {start-option.end()}")
        if start - option.end() > 6:
            if verbose:
                print("Found one SD valid motif")
            return True
    if verbose:
        print("No SD sequence found")
    return False


def predict_genes(sequence, start_regex, stop_regex, shine_regex,
                  min_gene_len, max_shine_dalgarno_distance, min_gap):
    """Predict most probable genes
    Based on teacher's algorithml written in French
    """
    #print(f"Studying a {len(sequence)} bases long sequence")
    predicted_genes = []

    start = 0
    while len(sequence) - start >= min_gap:
        start = find_start(start_regex, sequence, start, len(sequence))
        #print(f"starting position {start}")
        if start is None:
            break
        stop = find_stop(stop_regex, sequence, start)
        #print(f"found stop position {stop}")
        if stop is None:
            start += 1
            continue
        #print(f"current length {stop - start + 1} vs {min_gene_len}")
        if stop - start + 1 <= min_gene_len:
            # I would seek another stop but teacher's algo drop this start
            start += 1
            continue
        #sd_present = has_shine_dalgarno(shine_regex, sequence, start, max_shine_dalgarno_distance)
        #print(f"detected sd sequence: {sd_present}")
        if not has_shine_dalgarno(shine_regex, sequence, start, max_shine_dalgarno_distance):
            start += 1
            continue
        last_base = stop + 2 + 1  # +2 is 3rd codon letter, +1 for 1-based count
        predicted_genes.append([start+1, last_base])
        #print(f"saved gene positions {predicted_genes[-1]}")
        start = last_base + min_gap
        #print(f"start for next iteration: {start}")
    return predicted_genes


def write_genes_pos(predicted_genes_file, probable_genes):
    """Write list of gene positions
    """
    try:
        with open(predicted_genes_file, "wt") as my_file_out:
            predict_genes_writer = csv.writer(my_file_out, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
    except IOError:
        sys.exit("Error cannot open {}".format(predicted_genes_file))


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_genes(fasta_file, sequence, probable_genes, sequence_rc, probable_genes_comp):
    """Write gene sequence in fasta format
    """
    try:
        with open(fasta_file, "wt") as fasta:
            for i,gene_pos in enumerate(probable_genes):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                    i+1, os.linesep,
                    fill(sequence[gene_pos[0]-1:gene_pos[1]])))
            i = i+1
            for j,gene_pos in enumerate(probable_genes_comp):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                            i+1+j, os.linesep,
                            fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])))
    except IOError:
        sys.exit("Error cannot open {}".format(fasta_file))


def reverse_complement(kmer):
    """Get the reverse complement"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in kmer[::-1]])


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """

    def call_me_predict_genes(sequence):
        """Wrap predict_genes function with given main arguments.
        """
        min_length = args.min_gene_len
        max_sd_dist = args.max_shine_dalgarno_distance
        min_spacing = args.min_gap

        return predict_genes(sequence, start_regex, stop_regex, shine_regex,
                             min_length, max_sd_dist, min_spacing)

    # Gene detection over genome involves to consider a thymine instead of
    # an uracile that we would find on the expressed RNA
    #start_codons = ['TTG', 'CTG', 'ATT', 'ATG', 'GTG']
    #stop_codons = ['TAA', 'TAG', 'TGA']
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    # Shine AGGAGGUAA
    #AGGA ou GGAGG
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    # Arguments
    args = get_arguments()

    # Retrieve sequence
    sequence = read_fasta(args.genome_file)
    print(f"Input sequence length: {len(sequence)}")
    sequence.replace("U", "T")

    # Let us do magic in 5' to 3'
    print("Scanning forward sequence...")
    forward_genes = call_me_predict_genes(sequence)
    print("done")

    # Let's do the same in 3' to 5'
    print("Reverse complement generation")
    sequence_rc = reverse_complement(sequence)
    print("Scannning reverse sequence...")
    reverse_genes = call_me_predict_genes(sequence_rc)
    print("done")

    # Update reverse positions in 5'-3' values and merge them with forward positions
    print("Building final gene list")
    probable_genes = forward_genes[:]
    for end, start in reverse_genes[::-1]:
        end = len(sequence) - end
        start = len(sequence) - start
        probable_genes.append([start, end])
    print(f"{len(probable_genes):10d} probable genes found")
    probable_genes.sort()

    # Call to output functions
    print("Writting output files")
    write_genes_pos(args.predicted_genes_file, probable_genes)
    write_genes(args.fasta_file, sequence, forward_genes, sequence_rc, reverse_genes)


if __name__ == '__main__':
    main()
