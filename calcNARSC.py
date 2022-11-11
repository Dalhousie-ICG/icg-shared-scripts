#!/usr/bin/env python
import sys
from statistics import mean
from Bio import SeqIO
from Bio.Seq import Seq

# load fasta as file obj
my_fasta = open(sys.argv[1])

# state number of nitrogen
# atoms per amino acid
n_sc_counts = {
    'H': 2,
    'K': 1,
    'N': 1,
    'Q': 1,
    'R': 3,
    'W': 1,
}


# define function that calculates narsc per protein
# protein is a Seq() object,
# mapping a dict with atom counts per aa
def get_narsc_protein(protein, mapping=n_sc_counts):
    n_atom_count = 0
    for aa in mapping.keys():
        n_atom_count += protein.count(aa) * mapping.get(aa)
    narsc = n_atom_count / len(protein)
    return narsc


# test function
assert get_narsc_protein(protein=Seq("KNQWRAAAAA"), mapping=n_sc_counts) == 0.7

# parse my_fasta and get list of seq objs
seq_objs = [record.seq for record in SeqIO.parse(my_fasta, format='fasta')]
# get list of narscs and get the mean
print(mean(map(get_narsc_protein, seq_objs)))
