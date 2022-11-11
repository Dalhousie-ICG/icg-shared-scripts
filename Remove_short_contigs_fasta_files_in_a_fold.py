#!/usr/bin/python

'''
Remove contigs shorter than MINsize for all the fasta files in the folder.
Input an int (MINsize) for minimal length, and the suffix of fasta files (ignore if suffix = .fasta), i.e. file_name.fasta
Output:  file_name_MINsize.fna

Usage: python %Prog MIN_Size [suffix_of_fasta_files (default: fasta)]

'''

import os,sys
import glob
import argparse

def reading_infile (infile, separator, n):
# Read in file, separate headers with sequences and
# store these info in a dictionary

    with open(infile) as I:
        dic = {}
        I = I.read().split(separator)[1:]
        for read in I:
            entry = read.split('\n')
            key = entry[0]
            seq = ''.join(entry[1:-1])
            if len(seq) >= n:
                dic[key] = read
        return dic

		
def writing(infile, n, dic):
# Write the resulted sequences into output file: file_name_MINsize.fna
 
    outfile = infile.replace('.fasta','_%s.fna'% str(n))
    with open(outfile, 'w') as out:
        for key, value in dic.items():
            line = ">%s" % (value)
            out.write(line)


if __name__ == '__main__':
    minsize = int(sys.argv[1])
    if sys.argv[2]:
        files = glob.glob(".%s" % sys.argv[2])
    else:
        files = glob.glob(".fasta")


    for file in files:
        print("Now check file: %s" % file)
        contig_dic = reading_infile(file, '>', minsize)
        out = writing(file, minsize, contig_dic)
			


