#!/usr/bin/env python3

import argparse
from Bio import SeqIO

# E180_S6_vp1_cat.fasta
"""

The program remove the description and keeps the sequence identifyer of a fasta file
record.description = "" remove the description
"""

parser = argparse.ArgumentParser(prog = 'remove_desc_arg.py [Options]')
parser.add_argument('--infasta', type=str,help= 'paste path to input fasta file', required=True)
parser.add_argument('--outfasta', type=str,help= 'paste  file to write corrected fasta file',required=True)
args = parser.parse_args()

input_fasta = args.infasta
out_fasta = args.outfasta


def remove_fasta_description(input_file, output_file):
    """
    Removes the description from each sequence in a FASTA file.

    Args:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file.
    """
    with open(input_file, 'r') as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    for record in records:
        record.description = ""  # Remove the description

    with open(output_file, 'w') as handle:
        SeqIO.write(records, handle, "fasta")

# Example usage:
#infasta = "E180_S6_vp1_cat.fasta"
#outfasta = "E180_S6_vp1_cor.fasta"

def main():
    remove_fasta_description(input_fasta, out_fasta)
if __name__ == "__main__":
    main()


