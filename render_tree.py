#!/usr/bin/env python

from ete3 import Tree, TreeStyle, TextFace, NodeStyle, AttrFace, faces
from ete3 import SVG_COLORS
import csv
import argparse
import re
import sys


# make command line interface
parser = argparse.ArgumentParser(
    description="""Root, color and render a tree in PNG, SVG or PDF.
                Roots with outgroup defined in the mapping file"""
)
parser.add_argument(
    "-t", 
    dest="tree_file",
    metavar="TREE_FILE",
    # required=True, 
    help="Tree in Newick format"
)
parser.add_argument(
    "-m", 
    dest="mapping_file", 
    metavar="MAPPING_FILE",
    # required=True, 
    help="Tab separated mapping file. "
         "TaxonPattern \\t Clade \\t Color. "
)
parser.add_argument(
    "--title", 
    dest="title", 
    metavar="TITLE",
    # required=True, 
    help="Desired title of your tree figure "
)
parser.add_argument(
    "-o",
    dest="outfile",
    metavar="OUTFILE",
    # required=True,
    help="Name of the outfile. "
         "Extension determines format. "
         "For example file.pdf will output PDF format"
)
parser.add_argument(
    "--colors",
    dest="colors",
    required=False,
    action="store_true",
    help="Print all colors available in ETE3"
)
parser.add_argument(
    "--no-root",
    dest="noroot",
    required=False,
    action="store_true",
    help="Do not root the tree"
)
args = parser.parse_args()


# print all available colors if requested
if args.colors:
    print("\nColors you can use in your mapping file:\n")
    print( ', '.join(sorted( SVG_COLORS )) )
    print("\nGoogle SVG colors to check their actual colors")
    sys.exit()


# load tree
t = Tree(args.tree_file)
leaf_names = t.get_leaf_names()


# load mapping file
with open(args.mapping_file, "r") as f:
    map = {}
    for taxon, clade, color in csv.reader(f, delimiter="\t"):
        # find taxa in tree that match TaxonPattern in map file
        pattern = re.compile( taxon )
        for tree_taxon in filter(pattern.search, leaf_names):
            # assign tree_taxon a clade and a color 
            map[tree_taxon] = [clade, color]


# root tree
if not args.noroot:
    tree_outgroup_taxa = [taxon for taxon in map.keys() if map[taxon][0] == "Outgroup"]

    # in case the provided tree has the root already within the outgroup,
    # root outside the outgroup first
    tree_ingroup_taxa =  [taxon for taxon in map.keys() if map[taxon][0] != "Outgroup"]
    random_ingroup_taxon = tree_ingroup_taxa[0]
    t.set_outgroup(random_ingroup_taxon)

    # then root with the outgroup
    outgroup_lca = t.get_common_ancestor(tree_outgroup_taxa)
    t.set_outgroup(outgroup_lca)


# set TreeStyle
ts = TreeStyle()
ts.scale = 200  # 400 pixels per branch length unit
# ts.branch_vertical_margin = 5  # 5 pixels between adjacent branches
ts.show_leaf_name = False
ts.show_branch_support = True
ts.scale_length = 0.3
# # title = re.sub("_allTaxa_.*", "\n", args.tree_file)
# title = args.title
# ts.title.add_face(TextFace(title, fsize=20), column=0)


# set NodeStyle
ns = NodeStyle()
ns["size"] = 0 # suppress node icon
ns["hz_line_width"] = 2
ns["vt_line_width"] = 2
for node in t.traverse():
    node.set_style(node_style=ns)


# Order tree
t.ladderize(direction=1)


# Set leaf font
def leaf_font(node):
    """
    Set font for leaf names
    and replace underscores with spaces
    """
    if node.is_leaf():

        # extract color from loaded mapping file
        ## if node.name not present in mapping file, give it color black
        color = map.get(node.name, ['undefined_clade', 'black'])[1]

        # remove underscores from leaf name
        # node.name = node.name.replace("_", " ")

        # color leaf name
        leaf_face = AttrFace(attr="name", fgcolor=color, fsize=10)
        faces.add_face_to_node(
            face=leaf_face, node=node, column=1, position="branch-right"
        )


# render tree
t.render(args.outfile, tree_style=ts, layout=leaf_font, dpi=300, h=280, units="mm")
