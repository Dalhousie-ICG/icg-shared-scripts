#!/usr/bin/env python
"""
Copyright 2023 Jason Shao & Joran Martijn.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.
If not, see <https://www.gnu.org/licenses/>.

-----------------------------------------------------------------------------------------------------------------------

This program creates a simple unrooted tree from a newick string, makes BarChartFace for each leaf node in layout and
render said tree to a PNG image.

mandatory arguments: -t or --tree, newick tree
                     -n or --filename, path to alignment file
                     -f or --format, format of said alignment file (only tested for fasta)
optional arguments: -o or --output, desired name for output file, defaults to "tree".png
                    -s or --subsets, two subsets of amino acids to be highlighted, entered as a comma delimited string
                    -m or --frequency_type, type of frequency to display, "absolute" or "relative"
                    -g or --outgroup_reps, chosen outgroup(s) for rooting
                    -c or --show_chi2_score, whether to display chi-square score
"""
import argparse
import sys
import os
from Bio import SeqIO
from ete3 import Tree, faces, TreeStyle, BarChartFace, TextFace


# This class implements a holder for taxon attributes:
# name, seq, amino acid content, percentages, and the ability to calculate said attributes.
class Taxon:
    all_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    def __init__(self, name, seq):
        self.name = name
        self.seq = str(seq).replace("-", "")
        self.freqs = {aa: self.seq.count(aa) / len(self.seq) for aa in self.all_amino_acids}
        self.display_freqs = [self.freqs]
        self.display_max_value = 0.2
        self.chi_square_score = 0

    def set_all_relative_freq(self, avg_freq_dict):
        relative_freqs = {aa: self.freqs[aa] - avg_freq for aa, avg_freq in avg_freq_dict.items()}
        self.display_freqs = [relative_freqs]

    def set_subset_abs_freq(self, subsets):
        group1_freq = {}
        group2_freq = {}
        other_freq = {}

        for key, value in self.freqs.items():
            if key in subsets[0]:
                group1_freq[key] = value
            elif key in subsets[1]:
                group2_freq[key] = value
            else:
                other_freq[key] = value

        self.display_freqs = [group1_freq, group2_freq, other_freq]

    def set_subset_relative_freq(self, subsets, avg_freq_dict):
        group1_rel_freq = {}
        group2_rel_freq = {}
        other_rel_freq = {}

        for aa, avg_freq in avg_freq_dict.items():
            if aa in subsets[0]:
                group1_rel_freq[aa] = self.display_freqs[0][aa] - avg_freq
            elif aa in subsets[1]:
                group2_rel_freq[aa] = self.display_freqs[1][aa] - avg_freq
            else:
                other_rel_freq[aa] = self.display_freqs[2][aa] - avg_freq

        self.display_freqs = [group1_rel_freq, group2_rel_freq, other_rel_freq]

    def calculate_chi_square(self, align_freqs):
        chi_square_score = 0
        taxon_seq_len = len(self.seq)

        for aa in self.all_amino_acids:
            expected_count = align_freqs[aa] * taxon_seq_len
            observed_count = self.freqs[aa] * taxon_seq_len
            chi_square_score += (observed_count - expected_count) ** 2 / expected_count

        self.chi_square_score = round(chi_square_score, 1)


# check if subset is in the right format and if it contains valid amino acids
def validate_subsets(subset):
    subsets = [aa_group.upper() for aa_group in subset.split(",")]

    if len(subsets) != 2:
        raise argparse.ArgumentTypeError(f"Subsets does not have exactly 2 groupings: {subset}")

    if any(char not in all_amino_acids for aa_group in subsets for char in aa_group):
        raise argparse.ArgumentTypeError(f"At least one subset contains invalid amino acid(s)")

    return subsets


# check if frequency_type is valid
def validate_frequency(freq_type):

    if freq_type not in ["absolute", "relative"]:
        raise argparse.ArgumentTypeError(
            "Invalid tag for frequency type, only \"absolute\" or \"relative\" are allowed!"
        )

    return freq_type


# check if the provided outgroup is valid
def validate_outgroup(outgroup_reps, leaves):
    outgroup_reps = outgroup_reps.split(",")

    if any(rep not in leaves for rep in outgroup_reps):
        raise argparse.ArgumentTypeError("Invalid outgroup, make sure to check spelling!")

    return outgroup_reps


"""
the above functions are for argument validation
"""


# layout function pre-defined by ete3
def layout(node):

    # only consider leaf nodes
    if not node.is_leaf():
        return

    # retrieve Taxon from memory
    taxon = taxa_dict[node.name]
    dict_list = taxon.display_freqs

    i = 1
    for freq_dict in dict_list:
        face = get_barchart_face(freq_dict, taxon.display_max_value)

        # all faces have empty labels except for the bottom row
        root_node = node.get_tree_root()
        if node.name == root_node.get_leaf_names()[-1]:
            face.labels = list(freq_dict.keys())

        # ensure a healthy amount of space between the tree and the faces
        if i == 1:
            face.margin_left = 50
        face.margin_right = 10  # same for between faces

        # all faces have no scale except for the right-most column
        if len(freq_dict) > 1 and i != len(dict_list):
            face.scale_fsize = 1

        faces.add_face_to_node(face=face, node=node, column=i, position="aligned")
        i += 1

    # display chi-square scores if specified
    if taxon.chi_square_score != 0:
        text_face = TextFace(taxon.chi_square_score)
        text_face.margin_left = 50
        faces.add_face_to_node(face=text_face, node=node, column=i, position="aligned")


def get_barchart_face(freq_dict, max_value):
    face = BarChartFace(
        values=[abs(x) for x in freq_dict.values()],
        labels=[" " for x in freq_dict.keys()],
        label_fsize=9,  # this value dictates scaling if bar widths are uniform
        colors=["blue" if f > 0 else "red" for f in freq_dict.values()],
        width=40,  # when below a certain threshold, all the bar widths are scaled to be uniform
        height=50,
        max_value=max_value,
    )

    return face


# if user specified outgroup taxa in the flags then root accordingly
def root(tree, outgroup_reps):

    # check if the outgroup is only one taxon
    if len(outgroup_reps) == 1:
        tree.set_outgroup(outgroup_reps[0])  # just make that one taxon the outgroup

        return tree

    # if outgroup has multiple taxa then get their common ancestor
    common_ancestor = tree.get_common_ancestor(*outgroup_reps)

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


def main(args):

    global taxa_dict

    tree = Tree(args.tree, format=1)
    leaves = tree.get_leaf_names()
    outfile = args.output + ".png"

    #  check if outgroup(s) is specified
    if args.outgroup_reps is not None:

        try:
            outgroup_reps = validate_outgroup(args.outgroup_reps, leaves)
        except argparse.ArgumentTypeError as e:
            print(e)
            sys.exit()

        tree = root(tree, outgroup_reps)

    frequency_type = args.frequency_type
    subsets = args.subsets  # a list assumed to have two items
    show_chi2_score = args.show_chi2_score

    taxa_dict = {}  # dictionary to store taxa
    all_seq = ""  # string to hold all the sequences in the alignment

    # read in fasta to parse for Taxon objects
    for seq_record in SeqIO.parse(args.file, args.format):
        new_taxon = Taxon(seq_record.id, seq_record.seq)
        all_seq += str(seq_record.seq).replace("-", "")
        taxa_dict[seq_record.id] = new_taxon  # get taxa dict

    # calculate relative frequencies if specified
    if frequency_type == "relative" or show_chi2_score is True:
        # calculate average frequency
        avg_freq_dict = {aa: all_seq.count(aa) / len(all_seq) for aa in all_amino_acids}

    # determine which frequencies to display, as there are a total of 4 scenarios
    for taxon in taxa_dict.values():

        if subsets is None:
            # scenario 1: no subsets, absolute frequencies
            # not show as it is the default behaviour

            # scenario 2: no subsets, relative frequencies
            if frequency_type == "relative":
                taxon.set_all_relative_freq(avg_freq_dict)
                taxon.display_max_value = 0.05

        else:
            # scenario 3: subsets, absolute frequencies
            taxon.set_subset_abs_freq(subsets)

            # scenario 4: subsets, relative frequencies
            if frequency_type == "relative":
                taxon.set_subset_relative_freq(subsets, avg_freq_dict)
                taxon.display_max_value = 0.05

        if show_chi2_score is True:
            taxon.calculate_chi_square(avg_freq_dict)

    # tree styling
    tree.ladderize()
    tree_style = TreeStyle()
    tree_style.show_scale = False  # do not show scale

    # render tree
    tree.render(
        file_name=outfile,
        units="px", h=200 * len(leaves),
        tree_style=tree_style,
        layout=layout
    )

    print("\nFinished.\n")

    if os.path.exists(outfile):
        print(f"Tree is successfully created under {outfile}")
    else:
        print("Something went wrong, tree not generated.")


# Guard against undesired invocation upon import
if __name__ == "__main__":
    all_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    # specify options, disable for debugging
    parser = argparse.ArgumentParser(description="Tree making")

    parser.add_argument("-t", "--tree", required=True)
    parser.add_argument("-n", "--file", required=True)
    parser.add_argument("-f", "--format", required=True)
    parser.add_argument("-o", "--output", type=str, default="tree")
    parser.add_argument("-s", "--subsets", type=validate_subsets, nargs="?")
    parser.add_argument("-m", "--frequency_type", type=validate_frequency, default="absolute")
    parser.add_argument("-g", "--outgroup_reps", type=str, nargs="?")
    parser.add_argument("-c", "--show_chi2_score", type=bool, default=False)

    taxa_dict = None

    arguments = parser.parse_args()
    main(arguments)
