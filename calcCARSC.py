#!/usr/bin/env python
import sys
from statistics import mean
from Bio import SeqIO
from Bio.Seq import Seq

# load fasta as file obj
my_fasta = open(sys.argv[1])

# state number of nitrogen
# atoms per amino acid
c_sc_counts = {
    'R': 4,
    'H': 4,
    'K': 4,
    'D': 2,
    'E': 3,
    'S': 1,
    'T': 2,
    'N': 2,
    'Q': 3,
    'C': 1,
    'U': 1,
    'G': 0,
    'P': 3,
    'A': 1,
    'V': 3,
    'I': 4,
    'L': 4,
    'M': 3,
    'F': 7,
    'T': 7,
    'W': 9
}


# define function that calculates carsc per protein
# protein is a Seq() object,
# mapping a dict with atom counts per aa
def get_carsc_protein(protein, mapping=c_sc_counts):
    c_atom_count = 0
    for aa in mapping.keys():
        c_atom_count += protein.count(aa) * mapping.get(aa)
    carsc = c_atom_count / len(protein)
    return carsc


# test function
assert get_carsc_protein(protein=Seq("KNQWRAAAAA"), mapping=c_sc_counts) == 2.7

# parse my_fasta and get list of seq objs
seq_objs = [record.seq for record in SeqIO.parse(my_fasta, format='fasta')]
# get list of narscs and get the mean
print(mean(map(get_carsc_protein, seq_objs)))
