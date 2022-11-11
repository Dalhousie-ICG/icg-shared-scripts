#!/usr/bin/env python
"""
Author: Jimmy Saw
Date: 2014-02-21
Last updated: 2014-06-07 (Jimmy)
Last updated: 2022-05-17 (Joran Martijn)

Description: This script can draw matching regions between assemblies based on MUMMER alignments.
             You can also have MEGAN/Blastp taxonomic path file included, along with prodigal ORFs
             to show the location of ORFs that will be colored by taxonomic origin. If tax file is not
             given, it will not be colored.

Usage example:
    compare_assemblies.py -r -r flye_ergo_contigs.fasta -q canu_ergo_contigs.fasta -m flyeREF_vs_canuQRY.nucmer.coords -l 10000 -n 1000 -i 99.5

To get MUMMER coordinates, run:

$ nucmer <reference.fasta> <query.fasta> -p <nucmer_ref_vs_quer>
$ show-coords -r -c -l -T <nucmer_ref_vs_quer.delta> > <nucmer_ref_vs_quer.coords>

Make sure you don't mix up reference and query when running this script!!

"""
import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Bio import SeqIO


### Read command line arguments ###

parser = argparse.ArgumentParser(description="This script plots similarity regions between two genomes after running nucmer. Enhanced version of MUMMER plot.")
parser.add_argument("-r", "--reference_fasta", required=True, help="Reference contigs. This is the FIRST file given when running nucmer.")
parser.add_argument("-q", "--query_fasta", required=True, help="Subject contigss. This is the SECOND file given when running nucmer.")
parser.add_argument("-m", "--mummer_file", required=True, help="nucmer .coord file")
parser.add_argument("-l", "--min_aln_length", required=True, help="Minimum alignment length to show")
parser.add_argument("-i", "--min_aln_id", required=True, help="Minimum alignment similarity to show")
parser.add_argument("-n", "--min_contig_length", required=True, help="Contigs below this length will be colored gray")

args = parser.parse_args()
g1seqs = SeqIO.parse(args.reference_fasta, "fasta") #generator object FastaIterator
g2seqs = SeqIO.parse(args.query_fasta, "fasta") #generator object FastaIterator
mmfile = args.mummer_file # canu_vs_flye.nucmer.coords
minalnlen = int(args.min_aln_length) 
minalnid  = float(args.min_aln_id)
minctglen = int(args.min_contig_length) 


### Parse both assemblies ###

g1seqdict = {}
g2seqdict = {}

for seqrec in g1seqs:
    g1seqdict[seqrec.id] = seqrec
# seqrec is the SeqRecord() object
# seqrec.id is a simple string containing the id of the sequence

for seqrec in g2seqs:
    g2seqdict[seqrec.id] = seqrec

sorted_g1 = []
g1assembly_size = 0

for k, v in g1seqdict.items():
    """
    Make tuple of contig name and lengths
    """
    sorted_g1.append((k, len(v.seq)))
    g1assembly_size += len(v.seq)
# k = seqrec.id, v = seqrec object, v.seq = the sequence itself, len(v.seq) the length of that sequence
# sorted_g1 is a list of tuples [(seqrec.id, seq length),()]

sorted_g2 = []
g2assembly_size = 0

for k, v in g2seqdict.items():
    """
    Make tuple of contig name and lengths
    """
    sorted_g2.append((k, len(v.seq)))
    g2assembly_size += len(v.seq)

largest_total = max([g1assembly_size, g2assembly_size])
# largest_total is the length of the longest total assembly size between the two compared
# this is useful for plotting (see below)

sorted_g1.sort(key = lambda x: x[1], reverse=True)
sorted_g2.sort(key = lambda x: x[1], reverse=True)
# sorted_g1 and sorted_g2 are the list of tuples sorted by sequence length (largest to smallest)

# build dict of g1 assembly contig names and padsizes
padding_g1 = {}
g1count = 0
g1_padsize = 0
# len(sorted_g1) returns the number of contigs of assembly g1
while g1count < len(sorted_g1):
    """
    Make padding sizes for sorted contigs (largest -> smallest) so that contigs
    can be placed along x axis, and MUMMER result file can still be read to draw
    rectangles after padding. 
    """
    contig_name   = sorted_g1[g1count][0]
    contig_length = sorted_g1[g1count][1]
    ## sorted_g1[0][0] is the largest contig name
    ## sorted_g1[1][0] is the second largest contig name, etc
    padding_g1[contig_name] = g1_padsize  # put padding size into a dictionary of contig names
    ## [largest contig name]         = padsize 0
    ## [second largest contig name]  = padsize 2 136 038
    g1_padsize += contig_length
    g1count += 1

# build dict of g2 assembly contig names and padsizes
padding_g2 = {}
g2count = 0
g2_padsize = 0
while g2count < len(sorted_g2):
    contig_name   = sorted_g2[g2count][0]
    contig_length = sorted_g2[g2count][1]
    padding_g2[contig_name] = g2_padsize
    g2_padsize += contig_length
    g2count += 1


### Draw the contigs/scaffolds ###

fig = plt.figure(num=1, figsize=(20, 6))
# plt.figure(num=int) is unique identifier for the figure
# plt.figure(figsize=(float,float)) is width and height in inches
ax1 = fig.add_subplot(111)
# .add_subplot() is a method of plt.figure
# .add_subplot(111) is the same as .add_subplot(nrows=1,ncols=1,index=1)
# which describes the position of the subplot
# it will take the <index> position on a grid with <nrows> rows and <ncols> cols

for index, a in enumerate(sorted_g1):

    # get start and stop paddings from padding_g1 dict
    ## index = list index of tuple
    ## a = the tuple, i.e. ('contig_name', 'contig_length')
    contig_name   = a[0]
    contig_length = a[1]
    padded_start = padding_g1[contig_name]
    padded_stop  = padding_g1[contig_name] + contig_length

    # set y coordinates and color for contigs
    x = [padded_start, padded_stop]
    y = [2.05, 2.05]
    ## if index % 2.0 == 0, then index is even
    ## first contig y = 5.05, second contig y = 4.95, then 5.05, etc
    if index % 2.0 == 0:  # set alternating y coords to see contig boundaries
        y = [1.95, 1.95]
    ## contigs bigger than 1kb are black
    if contig_length >= minctglen:  
        ax1.plot(x, y, linewidth=4.0, color='#000000') # #000000 = black
    ## contigs smaller than 1kb are grey
    else:
        ax1.plot(x, y, linewidth=4.0, color='#AAAAAA') # #AAAAAA = grey
#     #ax1.annotate(a[0], xy=(padded_start, 5), fontsize=8, rotation=90)

for index, a in enumerate(sorted_g2):

    contig_name   = a[0]
    contig_length = a[1]
    padded_start = padding_g2[contig_name]
    padded_stop = padding_g2[contig_name] + contig_length

    x = [padded_start, padded_stop]
    y = [9.05, 9.05]
    if index % 2.0 == 0:
        y = [8.95, 8.95]
    if contig_length >= minctglen:
        ax1.plot(x, y, linewidth=4.0, color='#000000') # black
    else:
        ax1.plot(x, y, linewidth=4.0, color='#AAAAAA') # grey


### Draw the matching regions ###

# these are used to calculate total alignment 
# length between the two assemblies,
# which will reported in the plot
g1totalmatch = 0
g2totalmatch = 0


# parse nucmer coords file
## "r" open for reading
mummer = open(mmfile, "r")
mfl = mummer.readlines()
# mfl = a list of lines, where each line is a match
# mfl[4:] is the same list with the first 4 lines skipped
for line in mfl[4:]:

    # break line into list of fields
    c = line.strip().split('\t')

    ident = float(c[6]) # alignment similarity
    g1matchlen = int(c[4]) # alignment length from g1 perspective
    g2matchlen = int(c[5]) # alignment length from g2 perspective

    # skip alignment if identity is < 90%
    if ident < minalnid or g1matchlen < minalnlen:
        continue
    
    g1id = c[11] # g1 contig_name
    g2id = c[12] # g2 contig_name
    
    g1start = int(c[0]) # alignment start for g1
    g1stop  = int(c[1]) # alignment end   for g1
    g2start = int(c[2]) # alignment start for g2
    g2stop  = int(c[3]) # alignment end   for g2

    g1totalmatch += g1matchlen # add all alignment lengths together
    g2totalmatch += g2matchlen

    x1padding = padding_g1[g1id] # retrieve x padding for g1 contig from padding_g1 dict
    x2padding = padding_g2[g2id] # retrieve x padding for g2 contig from padding_g2 dict
    # calculate padded start and stop coordinates
    padded_g1start = g1start + x1padding
    padded_g1stop  = g1stop  + x1padding
    padded_g2start = g2start + x2padding
    padded_g2stop  = g2stop  + x2padding

    # only if alignment length is over minimum length threshold
    # draw plot
    if g1matchlen >= minalnlen:
        x = [padded_g1start, padded_g2start, padded_g2stop, padded_g1stop]
        y = [2, 9, 9, 2]
        ax1.fill(x, y, color="grey", alpha=0.5)

### Annotate plots ###

# plt.axis(xmin, xmax, ymin, ymax)
plt.axis( [(-0.01 * largest_total), largest_total + (0.01 * largest_total), 0, 12] )
plt.xlabel('Cumulative contig/scaffold length')
plt.annotate('REF = ' + os.path.basename(args.reference_fasta), xy=(g1assembly_size / 2.0, 1))
plt.annotate('QRY = ' + os.path.basename(args.query_fasta),     xy=(g1assembly_size / 2.0, 10))
plt.annotate('Reference assembly size = ' + str(g1assembly_size),  xy=(largest_total*0.01, 10.5), fontsize=10)
plt.annotate('Query assembly size = '     + str(g2assembly_size),  xy=(largest_total*0.01, 11.0), fontsize=10)
plt.annotate('Minimum alignment length = '+ str(minalnlen) + "bp", xy=(largest_total*0.01, 1.0),  fontsize=10)
plt.annotate('Minimum alignment id = '    + str(minalnid) + "%",   xy=(largest_total*0.01, 0.5),  fontsize=10)
#plt.annotate('Total alignment length for reference = ' + str(g1totalmatch), xy=(largest_total*0.01, 0.5),  fontsize=10)
#plt.annotate('Total alignment length for query = '     + str(g2totalmatch), xy=(largest_total*0.01, 10.5), fontsize=10)
ax1.set_yticklabels([])
ax1.set_yticks([])


outfile=os.path.splitext(mmfile)[0] + ".png"
plt.savefig(outfile, dpi=300)
plt.show()


#### Joran: superfluous code from Jimmy that may be useful  later on

    # determine fillcolor based on level of similarity
    # fillcolor = '#AAAAAA' # grey by default
    # if ident >= 99:
    #     fillcolor = "dimgrey"
    # elif ident >= 98:
    #     fillcolor = "grey"
    # elif ident >= 97:
    #     fillcolor = "darkgrey"
    # elif ident >= 96:
    #     fillcolor = "lightgrey"
    # elif ident >= 95:
    #     fillcolor = "whitesmoke"
    # # else:
    # #     fillcolor = DIFFCOLORS[5]

### Draw color legend ###

# pcx = largest_total-(largest_total * 0.1)
# pctlgndx = [pcx, pcx, pcx, pcx, pcx, pcx]
# pctlgndy = [4.00, 3.75, 3.5, 3.25, 3.0, 2.75]
# pctfc = ["dimgrey", "grey", "darkgrey", "lightgrey","whitesmoke"]
# pcttx = ['>=99', '>=98%', '>=97%', '>=96%', '>=95%']
# for a, b, c, d in zip(pctlgndx, pctlgndy, pctfc, pcttx):
#     #prect = Rectangle((a-42000, b), 40000, 0.25, facecolor=c, alpha=0.5)
#     prect = Rectangle((a-(largest_total*0.01), b), largest_total*0.01, 0.25, facecolor=c, alpha=0.5)
#     plt.gca().add_patch(prect)
#     ax1.annotate(d, xy=(a+(a*0.003), b+(b*0.01)), fontsize=7)
# # plt.annotate("Matching colors", xy=(pcx-(pcx*0.1), 1.25), fontsize=8)
# parser.add_argument("-o1", "--orfs1", help="ORFs for reference genome")
# parser.add_argument("-o2", "--orfs2", help="ORFs for query genome")
# parser.add_argument("-t1", "--tax1", help="MEGAN taxonomy file for reference")
# parser.add_argument("-t2", "--tax2", help="MEGAN taxonomy file for query")

### Color scales (dark to light) ###

# set colors
# DIFFCOLORS = plt.cm.Dark2(np.linspace(0, 1, 6))

# DCOLS = {'Archaea':'green', 'Bacteria':'blue', 'Eukaryota':'red', 'Viruses':'orange', 'Unknown':'gray'}
#DIFFCOLORS = plt.cm.afmhot(np.linspace(0, 1, 7))
#DIFFCOLORS = plt.cm.gist_earth(np.linspace(0, 1, 7))
#DIFFCOLORS = plt.cm.gist_ncar(np.linspace(0, 1, 7))
#DIFFCOLORS = plt.cm.gist_stern(np.linspace(0, 1, 7))
#DIFFCOLORS = plt.cm.ocean(np.linspace(0, 1, 7))
#DIFFCOLORS = plt.cm.terrain(np.linspace(0, 1, 7))
#DIFFCOLORS = plt.cm.Set1(np.linspace(0, 1, 7))
#DIFFCOLORS = plt.cm.spectral(np.linspace(0, 1, 7))
#DIFFCOLORS = plt.cm.Dark2(np.linspace(0, 1, 7))
#DIFFCOLORS = plt.cm.hsv(np.linspace(0, 0.5, 6))
#DIFFCOLORS = plt.cm.gist_ncar(np.linspace(0.1, 0.8, 6))
#DIFFCOLORS = plt.cm.cubehelix(np.linspace(0.2, 1, 7))
#DIFFCOLORS = plt.cm.CMRmap(np.linspace(0.3, 1, 7))
#DIFFCOLORS = plt.cm.YlOrBr(np.linspace(1, 0, 7))
#DIFFCOLORS = plt.cm.gist_earth(np.linspace(0, 1, 7))

# idlgndy = [i - 1.75 for i in pctlgndy]
# domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses', 'No hits']
# domainCols = ['green', 'blue', 'red', 'orange', 'gray']

# for a, b, c, d in zip(pctlgndx, idlgndy, domainCols, domains):
#     idrect = Rectangle((a - (largest_total*0.01), b), largest_total*0.01, 0.25, facecolor=c, alpha=0.5)
#     plt.gca().add_patch(idrect)
#     plt.annotate(d, xy=(a+(a*0.003), b+(b*0.03)), fontsize=7)
# plt.annotate("ORF colors", xy=(pcx-(pcx*0.1), 1.5), fontsize=8)

# ### Draw the ORFs ###

# tax1dict = {}
# tax2dict = {}

# if args.tax1:
#     with open(args.tax1, "rU") as tf1:
#         tf1lines = tf1.readlines()
#         for line in tf1lines:
#             c = line.strip().split(";")
#             tax1dict[c[0]] = c[1:]

# if args.tax2:
#     with open(args.tax2, "rU") as tf2:
#         tf2lines = tf2.readlines()
#         for line in tf2lines:
#             c = line.strip().split(";")
#             tax2dict[c[0]] = c[1:]

# unknown = ['No hits', 'Not assigned']
# hasroot = ['root']

# if args.orfs1:
#     orfs1 = [orf for orf in SeqIO.parse(args.orfs1, "fasta")]
#     for orf in orfs1:
#         orf_start = int(orf.description.split(" ")[2])
#         orf_stop = int(orf.description.split(" ")[4])
#         frame = int(orf.description.split(" ")[6])
#         p = orf_regex.match(orf.id)
#         parent_ctg = p.group(1)
#         c = "gray"
#         if orf.id in tax1dict:
#             taxinfo = tax1dict[orf.id]
#             if taxinfo[1].strip() in hasroot:
#                 if taxinfo[3].strip() == "cellular organisms":
#                     domain = taxinfo[5].strip()
#                     if domain in DCOLS:
#                         c = DCOLS[domain]
#                 elif taxinfo[3].strip() in DCOLS:
#                     c = DCOLS[taxinfo[3].strip()]
#         if parent_ctg in padding_g1:
#             padded_start = padding_g1[parent_ctg] + orf_start
#             padded_stop = padding_g1[parent_ctg] + orf_stop
#             if frame == 1:
#                 ax1.plot([padded_start, padded_stop], [5.02, 5.02], color=c, lw=2, alpha=0.8)
#             else:
#                 ax1.plot([padded_start, padded_stop], [4.98, 4.98], color=c, lw=2, alpha=0.8)

# if args.orfs2:
#     orfs2 = [orf for orf in SeqIO.parse(args.orfs2, "fasta")]
#     for orf in orfs2:
#         orf_start = int(orf.description.split(" ")[2])
#         orf_stop = int(orf.description.split(" ")[4])
#         frame = int(orf.description.split(" ")[6])
#         p = orf_regex.match(orf.id)
#         parent_ctg = p.group(1)
#         c = "gray"
#         if orf.id in tax2dict:
#             taxinfo = tax2dict[orf.id]
#             if taxinfo[1].strip() in hasroot:
#                 if taxinfo[3].strip() == "cellular organisms":
#                     domain = taxinfo[5].strip()
#                     if domain in DCOLS:
#                         c = DCOLS[domain]
#                 elif taxinfo[3].strip() in DCOLS:
#                     c = DCOLS[taxinfo[3].strip()]
#         if parent_ctg in padding_g2:
#             padded_start = padding_g2[parent_ctg] + orf_start
#             padded_stop = padding_g2[parent_ctg] + orf_stop
#             if frame == 1:
#                 ax1.plot([padded_start, padded_stop], [9.02, 9.02], color=c, lw=2, alpha=0.8)
#             else:
#                 ax1.plot([padded_start, padded_stop], [8.98, 8.98], color=c, lw=2, alpha=0.8)
