#!/usr/bin/env python

'''
Install required packages with
$ conda create -n <env_name> pyhmmer biopython ete3
'''

"""
Copyright 2023 Jason Shao & Joran Martijn.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.
If not, see <https://www.gnu.org/licenses/>.
"""

from pyhmmer.hmmer import hmmscan
from pyhmmer.easel import MSAFile
from pyhmmer.plan7 import HMMFile
from ete3 import Tree, SeqMotifFace, AttrFace, TreeStyle, NodeStyle, faces
from Bio import SeqIO
import itertools
import re
import argparse

# command line options
parser = argparse.ArgumentParser(
            description=
        '''
        Plot domain composition for each protein in a phylogenetic tree,
        next to that tree to a PNG image file.

        The script takes in an alignment in FASTA format,
        a tree in NEWICK format, a database of HMM profiles in .h3m format,
        and optionally a pair of outgroup representatives to root the tree.

        The script uses pyhmmer to perform hmmscan, hence HMMER is technically
        not required.
        ''')

parser.add_argument(
        "-a", "--aln",
        dest='aln_file',
        type=str, required=True,
        help="Input alignment file")

parser.add_argument(
        "-t", "--tree",
        dest='tree_file',
        type=str, required=True,
        help="Input tree file")

parser.add_argument(
        "-d", "--hmm_db",
        dest='hmm_db',
        type=str, required=True,
        help="Input HMM profile database")

parser.add_argument(
        "-r", "--outgroup_reps", 
        type=str, required=False,
        nargs="?",
        help="outgroup_taxon1,outgroup_taxon2")

parser.add_argument(
        "-l", "--highlight", 
        type=str, required=False,
        nargs="?",
        help="Regular expression for taxa to highlight")

parser.add_argument(
        "-u", "--ungapped_seqs", 
        action='store_true', required=False,
        help="Invoke to see domains on ungapped sequences")

parser.add_argument(
        "-o", "--out", 
        type=str, required=False,
        default="tree.png",
        help="Output image in PNG format")

args = parser.parse_args()

def transform_coordinate(gapped_seq, ungapped_coordinate):
    i = 0 # counter for all positions including '-'
    j = 0 # counter for all positions excluding '-'
    for c in gapped_seq:
        i += 1
        if c != '-':
            j += 1
            if j == ungapped_coordinate:
                return i
    # if i is not found
    return None

# if user specified outgroup taxa in the flags then root accordingly
def root(tree, outgroup_reps):
    # check if the outgroup is only one taxon
    if len(outgroup_reps) == 1:
        tree.set_outgroup(outgroup_reps[0])  # just make that one taxon the outgroup
        return tree

    # if outgroup has multiple taxa then get their common ancestor
    common_ancestor = tree.get_common_ancestor(*outgroup_reps)
    # print(outgroup_reps)
    # print(common_ancestor)

    # if common ancestor is not the root then re-root with it as outgroup
    if not common_ancestor.is_root():
        tree.set_outgroup(common_ancestor)
        return tree

    # this is a workaround as you cannot re-root to the current "root" with ete3
    # therefore the user needs to provide an in-group taxon
    print("\nOutgroup common ancestor is the root node, an in-group taxon is required for re-rooting")
    ingroup_taxon = input("\nEnter an ingroup taxon: ")

    # reject non-leaf inputs
    while ingroup_taxon not in tree.get_leaf_names():
        print("You entered a taxon that is not a part of the tree, check spelling!")
        ingroup_taxon = input("\nEnter an ingroup taxon: ")

    # reject outgroup inputs
    while ingroup_taxon in outgroup_reps:
        print("You entered an outgroup taxon, check spelling!")
        ingroup_taxon = input("\nEnter an ingroup taxon: ")

    # setting outgroup with an in-group taxon ensures the real outgroup common ancestor is not the root node
    tree.set_outgroup(ingroup_taxon)
    common_ancestor = tree.get_common_ancestor(*outgroup_reps)
    tree.set_outgroup(common_ancestor)

    return tree

def leaf_font(node):
    # only consider leaf nodes
    if not node.is_leaf():
        return

    # default color is black
    color = "black"

    # highlight taxa that match regular expression red
    if args.highlight:
        pattern = re.compile('.*' + args.highlight)
        if re.match(pattern, node.name):
            color = "red"

    # color leaf name
    leaf_face = AttrFace(attr="name", fgcolor=color, fsize=10)
    faces.add_face_to_node(
        face=leaf_face, node=node, column=1, position="branch-right"
    )

if __name__ == '__main__':
    # store alignment in memory as SeqIO
    fasta = SeqIO.index(args.aln_file, 'fasta')

    # store alignment file in memory as MSAFile
    with MSAFile(args.aln_file, format="afa", digital=True) as msa_file:
        # .read() returns an MSA objects,
        # in this case a DigitalMSA object
        # msa.read().sequences returns a _DigitalMSASequences object
        # which is an iterable of DigitalSequence objects
        # to make it parseable by hmmscan(),
        # convert it to a list of DigitalSequence objects??
        seqs = list( msa_file.read().sequences )


    taxa2motifs = {}
    domains2colors = {}
    color_schemes = itertools.cycle( [ "rgradient:orange", 
                                       "rgradient:blue", 
                                       "rgradient:pink", 
                                       "rgradient:green", 
                                       "rgradient:purple" ] )

    # load hmm database into memory
    with HMMFile(args.hmm_db) as hmm_file:
        hmms = list(hmm_file)

    # store HMM lengths in dictionary
    hmm_lengths = { h.name.decode() : h.M for h in hmms }

    # execute hmmscan
    for query in hmmscan(seqs, hmms, cpus=0, E=0.001, domE=0.01):

        # skip query if it has no hits
        # satisfying the reporting E-value threshold
        if not query.reported:
            continue

        # get starts and ends of domains
        # and generate motifs for each domain
        motifs = []
        for hit in query.reported:
            for domain in hit.domains.reported:
                domain_name = hit.name.decode()

                # if its a newly identified domain
                # link a new colorscheme to it in memory
                if domain_name not in domains2colors.keys():
                    # grab color from memory
                    color = next( color_schemes )
                    domains2colors[domain_name] = color

                # retrieve aln seq from memory
                gapped_seq = str( fasta[ query.query_name.decode() ].seq )

                if not args.ungapped_seqs:
                    # transform ungapped_seq coordinate to gapped_seq coordinate
                    gap_from = transform_coordinate(gapped_seq, domain.env_from)
                    gap_to   = transform_coordinate(gapped_seq, domain.env_to)
                else:
                    gap_from = domain.env_from
                    gap_to   = domain.env_to

                # if N-terminus is missing,
                # mark domain with a red bar
                if domain.alignment.hmm_from > 30:
                    truncation_motif = [ gap_from + 1, gap_from + 3,
                                         "[]", 100, 10, "black", "black", None ]
                    motifs.append(truncation_motif)

                # if C-terminus is missing,
                # mark domain with a red bar
                target_length = hmm_lengths[domain_name]
                if target_length - domain.alignment.hmm_to > 30:
                    truncation_motif = [ gap_to - 3, gap_to - 1,
                                         "[]", 100, 10, "darkred", "darkred", None ]
                    motifs.append(truncation_motif)

                # generate motif using domain coordinates
                # and custom domain colors
                motif = [ gap_from, gap_to,
                          "[]", 100, 10, "black", domains2colors[domain_name],
                          f"arial|7|black|{domain_name}" ]
                motifs.append(motif)

        # link motifs to taxon (the query)
        taxa2motifs[ query.query_name.decode() ] = motifs

    # load tree
    tree = Tree(args.tree_file)

    # root tree
    if args.outgroup_reps:
        tree = root(tree, args.outgroup_reps.split(",") )

    # add faces to leaves
    for leaf in tree.get_leaves():
        # get unaligned sequence from aligned FASTA file
        if not args.ungapped_seqs:
            seq = str( fasta[leaf.name].seq )
        else:
            seq = str( fasta[leaf.name].seq.replace('-','') )

        # generate SeqMotifFace and attach to leaf node
        face = SeqMotifFace(seq=seq, seq_format="line", 
                            gap_format="blank", motifs=taxa2motifs[leaf.name])
        # face.opacity = 0.5
        leaf.add_face(face, 0, "aligned")

    # define tree style
    tree.ladderize(direction=1)
    tree_style = TreeStyle()
    tree_style.show_scale = False
    tree_style.show_leaf_name = False

    # define node style
    ns = NodeStyle()
    ns["size"] = 0
    ns["hz_line_width"] = 1
    ns["vt_line_width"] = 1
    for node in tree.traverse():
        node.set_style(node_style=ns)

    # render tree
    tree.render(args.out, h=50 * len(tree.get_leaf_names()), tree_style=tree_style, layout=leaf_font)
