#!/usr/bin/env python3
# to care about opened reading frame, split sequence into triplets and use itertools.takewhile to handle that but it can be done using loops as well:

import argparse
from Bio import SeqIO
from itertools import takewhile
import itertools 
parser = argparse.ArgumentParser(prog = 'dna2proteins_v2.py [Options]')
parser.add_argument('--fastapath', type=str,help= 'paste path to raw sample.fasta file', required=True)
parser.add_argument('--outpro', type=str,help= 'paste txt file to write ID protein fasta file',required=True)
args = parser.parse_args()
fastafile_path = args.fastapath
out_prot = args.outpro
#opening the fasta file
"""
def read_fasta(file_name):
    sequences = {}
    with open(file_name, 'r') as file:
        header = ""
        sequence = ""
        for line in file:
            if line.startswith('>'):
                if header:
                    sequences[header] = sequence
                header = line.strip()[1:]  # Remove the '>'
                sequence = ""
            else:
                sequence += line.strip()
        if header:  # Add the last sequence
            sequences[header] = sequence
    return sequences

for header, seq in sequences.items():
    print(header)
    print(seq)

"""

def read_fasta_biopython(file_name):
    sequences = []
    for record in SeqIO.parse(file_name, "fasta"):
        sequences.append(record)
    return sequences

#sequences = read_fasta_biopython("your_fasta_file.fasta")

#def print_seq(sequences):
    #for record in sequences:
        #print(record.id)
        #sequence = record.seq
        #return sequence
# codon dictionary

codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

#def translate_dna(sequence, codontable, stop_codons = ('TAA', 'TGA', 'TAG')):

def translate_dna(sequence, codontable, stop_codons = ('TAA', 'TGA', 'TAG')):       
    start = sequence.find('ATG')

    # Take sequence from the first start codon
    trimmed_sequence = sequence[start:]

    # Split it into triplets
    #codons = [trimmed_sequence[i:i+3] for i in range(0, len(trimmed_sequence), 3)]
    codons = [trimmed_sequence[i:i+3] for i in range(0, len(trimmed_sequence), 3)]
    #print(len(codons))
    #print(trimmed_sequence)
    #print(codons)

    # Take all codons until first stop codon
    coding_sequence  =  takewhile(lambda x: x not in stop_codons and len(x) == 3 , codons)

    # Translate and join into string
    protein_sequence = ''.join([codontable[codon] for codon in coding_sequence])

    # This line assumes there is always stop codon in the sequence (the _ is place of stop codon. 0 is 1st item in the parentheses)
    #return "{0}_".format(protein_sequence)
    return "{0}".format(protein_sequence)
# the codons = [trimmed_sequence[i:i+3] for i in range(len(trimmed_sequence)/3)]
# should actually look like codons = [trimmed_sequence[i:i+3] for i in range(0, len(trimmed_sequence), 3)]

"""
def add_header_see_fasta(fasta_file, header_list):
    with open(fasta_file, 'w') as f:
        for i, line in enumerate(lines):
            if line.startswith('>'):  # Check if the line is a header
                f.write(f">{header_list[i]}\n")
            else:
                f.write(line)
"""
#def add_sequence_and_header(file_path, header, sequence):


"""
def add_sequence_and_header(out_prot, header, sequence):
    with open(out_prot, 'w') as fi:
        #fi.write(f"{header}\n{sequence}\n")
        #f.write(header + '\n')
        #fi.writelines('\n'.join(header))
        fi.writelines('\n'.join(sequence))
"""
# combine the header list and sequence list
def add_header2sef(out_prot, header, sequence):
     jj = list(itertools.chain(*zip(header,sequence)))
     with open(out_prot, 'w') as fi:
         fi.writelines('\n'.join(jj))
   

"""
# Example usage
file_path = 'samp_protein.fasta'
header = '>'+ record.id
#sequence = 'ATCGTACG'

add_sequence_and_header(file_path, header, sequence)
"""

def main():
    sequences = read_fasta_biopython(fastafile_path)
    #sequence = print_seq(sequences)
    #protein_seq = translate_dna(sequences,codontable,stop_codons = ('TAA', 'TGA', 'TAG'))
    #print(sequences)
    #print(protein_seq)
    list_recid = [] # to store all record id
    list_protseq =[]
    for record in sequences:
        #print(record.id)
        #print(record.seq)
        protein_seq = translate_dna(record.seq,codontable,stop_codons = ('TAA', 'TGA', 'TAG'))
        #print(protein_seq)
        #print(record.id)
        #list_recid.append(record.id) 
        header = '>'+ record.id
        list_recid.append(header)
        list_protseq.append(protein_seq) 
        sequence = protein_seq
    #for header_prot, seqprot in zip(list_recid, list_protseq):
    prot_head  = list(zip(list_recid, list_protseq))
    #for tup in prot_head:
    #add_sequence_and_header(out_prot, tup[0], tup[1])
    #add_sequence_and_header(out_prot, list_recid,list_protseq)
    add_header2sef(out_prot, list_recid, list_protseq)
    #print(tup[0])
    #print(tup[1])
    print(list_protseq)
    #for i in list_protseq:
    #add_sequence_and_header(list_protseq)
    
    #print(list_recid)
    #print(prot_head)
    #fi.close()    
if __name__ == "__main__":
    main()

#fi.close()
