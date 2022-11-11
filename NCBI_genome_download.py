#!/usr/bin/python
__author__ = 'D.S.L, D.Z.'
__version__ = 'Nov_2018_v1.2'

'''
Usage:    python3 Genome_download.py <genome_IDs_to_recover.txt>

genome_IDs_to_recover.txt:  a txt file with genome accession numbers to be downloaded.
(genome accession without any version, see example below:
GCA_000002725
GCA_000002765
GCA_900500625)

You can use the following commands to get rid of version number:
cat genome_IDs.txt | cut -d "." -f1 > genome_IDs_to_recover.txt

This script will automatically download the newest version of NCBI assembly report from
ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt

'''


import sys, os
import subprocess 
from os import path 
from subprocess import PIPE, Popen



def my_run(cmd): 
    unfound_list=[]
    print("\tRun command: %s\n" %cmd)
    task = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True) 
    output = task.stdout.decode('UTF-8').split('\n')
    output_strip =[f.strip() for f in output if f.strip() != '']
    #print("\nLines in output_strip in my_run:\n")
    for line in output_strip: 
        print(line)  

def readingCatalog(infile):
    with open(infile) as In:
        FTP = {}
        for line in In:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                acc = line[0].split('.')[0]
                FTP[acc] = line[19]
        #print("Keys in FTP[0:3]: \n%s" % list(FTP.keys())[0:3])
        #print("Values in FTP[0:3]: \n%s" % list(FTP.values())[0:3])
        return FTP

def readgenomeIDs(infile):
    with open(infile) as I:
        return {line.strip() for line in I if not line == ''}

def perform(Catalog,IDs):
    for id in IDs:
        print('about to retrieve %s'%id)
        basename = Catalog[id].split('/')[-1]
        print('basename', basename)
        ftpcmd = 'wget %s/%s_genomic.fna.gz' % (Catalog[id],basename)
        #print('ftpcmd', ftpcmd)
        try:
            retrieve = Popen(ftpcmd,stderr=PIPE,
                             stdout=PIPE,shell=True)
            o,e = retrieve.communicate()
            os.system('gunzip %s_genomic.fna.gz'%basename)
            #print(o)
        except:
            print('ftpcmd'%ftpcmd)
    return


if __name__ == '__main__':
    report_file = 'assembly_summary_genbank.txt'
    if not path.exists(report_file):
        dw_file = "wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
        my_run(dw_file)
        
    Genomic_Catalog= readingCatalog(report_file)
    Genomes2Recover = readgenomeIDs(sys.argv[1])
    print("In total %s genomes to be downloaded." % len(Genomes2Recover))
    perform(Genomic_Catalog,Genomes2Recover)
