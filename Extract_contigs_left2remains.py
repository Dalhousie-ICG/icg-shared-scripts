#!/usr/bin/python

"""
Extract configs from fasta file.                                                                               
python3 %Prog Input.fasta wanted_ID_list.txt                                                    
"""                                                                                 
import sys,re                                                                          
import argparse
from Bio import SeqIO

def Main(fasta_file, wanted_IDs,suffix):
    basename = wanted_IDs.split(".txt")[0]
    # Output1: wanted contigs
    result_file1 = fasta_file.replace('.%s'%suffix,'_%s_extracted.%s'%(basename, suffix))
    # Output2: the remaining contigs
    result_file2 = fasta_file.replace('.%s'%suffix,'_remained.%s'%suffix)

    wanted = [line.strip() for line in open(wanted_IDs) if not line == '']  

    print("Extract %s contigs to %s" %(len(wanted), result_file1))
    print("Resulted remaining file: %s" %result_file2)
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    with open(result_file1, "w") as f , open(result_file2,"w") as O:
        for seq in fasta_sequences:
            if seq.id in wanted:

                SeqIO.write([seq], f, "fasta")
            else:
                SeqIO.write([seq], O, "fasta")				
    print("Finished")

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Extract configs from fasta file. \nResults: 1. fasta with required contigs, 2. fasta with remaining contigs.", 
                                     epilog='''Usage example:
                                           python3 Extract_contigs.py -i <fasta_file> -l <list_Ids>''',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i","--infile",help="Input fasta file", required=True)
    parser.add_argument("-l","--list",help="Input interesting sequence IDs, one per line, ending with '.txt'", required=True)


    args = parser.parse_args()
    fasta_file = args.infile
    list_file = args.list
    suffix = fasta_file.split(".")[-1]
    #print(suffix)
    Main(fasta_file, list_file, suffix)
