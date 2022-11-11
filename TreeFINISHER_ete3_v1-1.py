#!/bin/python3.8
# Yana Eglit
# v. 1.1, May 2022

import argparse
import sys
import csv
from os import path, getcwd
import glob
# import copy
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, AttrFace, faces, RectFace, StackedBarFace
from ete3.parser.newick import NewickError
import itertools
import re
import time  # for logfile

# intree = sys.argv[1]
# tax_to_fullname = sys.argv[2]
# taxgroups = sys.argv[3]
# colourlist = "ColourPalleteNew.tsv"

logfile = "TreeFINISHER.log"
tree_leaf_font_size = 10  # move into dictionary?

"""
Args to add: 
- default settings config?
- fill taxa not included in list with colours or not
    - for future: 'capture' list of taxa in clades that weren't included before (either specify only new, or all)
    (eg expanded Rhizaria sampling, request updated taxon list of Rhizaria)
- scaling: full, fit to page [select: A4, letter, custom in px]
- settings for showing full support as dots or numbers; showing numbers for support > X

To do:
- improve scalebar: text position, size, vertical bars (redo entirely?)
- finish logfile
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Processes input newick tree using Ete3 toolkit. Version 1.1. Unfinished and untested.
    """)
    parser.add_argument("-c", "--colour-list",
                        dest="colourlist",
                        type=str,
                        help="Tabbed file specifying colours; formatted as:\ngroup-name\t#colour-hex-code"
                        )
    parser.add_argument("-t", "--taxa-by-groups",
                        dest="taxgroups",
                        type=str,
                        help="File specifying which group each taxon belongs to.\
                         Expected formatting:\ntaxon\tgroup-name"
                        )
    parser.add_argument("-f", "--taxon-full-names",
                        dest="tax_to_fullname",
                        type=str,
                        help="File specifying desired expanded taxon names, formatted as:\ntaxon\tfull-name"
                        )
    parser.add_argument("-i", "--intree",
                        type=str,
                        # nargs="+",  # allows one or more input file, as a list
                        help="Input Newick tree. Use -b for multi bootstrap values.",
                        required=True
                        )
    parser.add_argument("-o", "--outfile",
                        type=str,
                        help="Specify output filename. If overwrite mode off, will return error if file exists.\
                         Default naming is {intree}{int}.svg",
                        required=False
                        )
    parser.add_argument("-w", "--overwrite",
                        action="store_true",
                        default=False,
                        help="Will overwrite as {intree}.svg or if --outfile specified, as {outfile}.",
                        required=False
                        )
    # parser.add_argument("-w"#skip, "--fill-in-taxa",
    #                     action="store_true",
    #                     default=False,
    #                     help="Colour taxa that are not in colour list but inside group",
    #                     # second argument: add options (same colour, red) (to indicate missing taxa in list)
    #                     required=False
    #                     )
    parser.add_argument("-s", "--coverage-stats",
                        dest="cov_stats",
                        type=str,  # do not use type=list!
                        nargs=2,
                        help="Adds coverage stats, number of genes and perc sites, from specified file,\
                        followed by total gene count. Eg. -s {stats.tsv} 254\n \
                        Input file must be formatted as follows: Taxon\tNumber_of_genes(int)\tPercent_sites(float)\
                        \nIf no genes plot desired, simply leave that column present but blank.",
                        required=False
                        )


args = parser.parse_args()

tax_to_fullname = args.tax_to_fullname
taxgroups = args.taxgroups
colourlist = args.colourlist
intree = args.intree


# if len(args.intree) == 1:  # for inputting a list of trees -- later
#     intree = args.intree[0]
#     multi_tree_mode = False
# else:
#     intree = args.intree
#     multi_tree_mode = True


def write_to_logfile(list_of_strings, lf=logfile):
    with open(lf, "a") as log:
        log.write("\n"+"\n".join(list_of_strings))


#####################
# Check input files #
#####################

def check_input_files(f):
    try:
        with open(f, "r") as fl:
            print(f"{f} exists;")
            first_line = fl.readline().rstrip()
    except FileNotFoundError:
        print(f"{f} not found!")
    return first_line


def check_cov_stats(f):  # later refactor to check whole columns instead and not mind extra data incl top header?
    """
    Checks formatting of first line of coverage stats file, and whether genes, sites, or both (or none)
    :param f: coverage stats file, tab-delimited
    :return: returns whether to include each of: gene_count, perc_sites
    """
    first_line = check_input_files(cov_stats_file).split("\t")
    try:
        if len(first_line[1]) == 0:  # checking second column, for gene count (should be int)
            print(f"Second column blank, skipping gene counts.")
            incl_gene_count = False
        else:
            int(first_line[1])
            incl_gene_count = True
    except ValueError:
        print(f"Is {f} formatted correctly? Expecting integer in second column. Skipping gene count.")
        incl_gene_count = False
    except IndexError:  # if columns 2-3 missing entirely
        print(f"Is {f} formatted correctly? Missing second and third columns. Proceeding without coverage stats.")
        incl_gene_count = False
        incl_perc_sites = False
        return [incl_gene_count, incl_perc_sites]  # this ends the function before the next try
    try:
        if len(first_line[2]) == 0:  # checking second column, for perc sites (should be float or int)
            print(f"Third column blank, skipping perc sites.")
            incl_perc_sites = False
            return [incl_gene_count, incl_perc_sites]  # exit function
        if 0 <= float(first_line[2]) <=100:
            incl_perc_sites = True  # if present and formatted correctly
        else:
            print(f"Is {f} formatted correctly? Expecting number between 0 and 100. Skipping perc sites.")
            incl_perc_sites = False
    except ValueError:
        print(f"Is {f} formatted correctly? Expecting number in second column. Skipping perc sites.")
        incl_perc_sites = False
    except IndexError:
        print(f"Missing third column. Proceeding without perc sites.")
        incl_perc_sites = False
    return [incl_gene_count, incl_perc_sites]


# check if coverage stats inputs correct, and that the tabbed file makes sense
if args.cov_stats:
    print(args.cov_stats)
    cov_stats_file = args.cov_stats[0]
    if args.cov_stats[1].isdigit():
        cov_stats_maxgenes = int(args.cov_stats[1])
    else:
        sys.exit(f"Incorrect max gene number input for --coverage-stats! Expected integer, got {args.cov_stats[1]}\n")
    incl_genes_sites = check_cov_stats(cov_stats_file)
    # if incl_genes_sites[0] == True:  # first column is gene count  # might not need this
    #     incl_cov_genes = True
    # if incl_genes_sites[1] == True:  # second column is perc sites
    #     incl_cov_sites = True
    if incl_genes_sites == [False, False]:
        args.cov_stats = None  # cancels cov_stats for the rest of the script  # probably a bad idea to do it this way


try:  # check if second column exists AND starts with #hex colour
    if check_input_files(colourlist).split("\t")[1].startswith("#"):
        pass
    else:
        print(f"Is the second column in {colourlist} formatted correctly as #COLOUR in hex?")
except IndexError:
    print(f"Does {colourlist} contain a second column for #COLOUR in hex?\n")
    print(f"First line is: {check_input_files(colourlist)}")


try:  # check if second column exists AND is a string
    if type(check_input_files(tax_to_fullname).split("\t")[1]) is str:
        pass
    else:
        print(f"Is the second column in {tax_to_fullname} formatted correctly as a full taxon name?")
except IndexError:
    print(f"Does {tax_to_fullname} contain a second column for full taxon names?\n")
    print(f"First line is: {check_input_files(tax_to_fullname)}")


try:  # check if second column exists AND is a string
    if type(check_input_files(tax_to_fullname).split("\t")[1]) is str:
        pass  # count number of taxon groups
        # check if each of them exists in colour list! (see function below)
    else:
        print(f"Is the second column in {taxgroups} formatted correctly ")
except IndexError:
    print(f"Does {taxgroups} contain a second column for full taxon names?\n")
    print(f"First line is: {check_input_files(taxgroups)}")

# if args.weirdboots:
#     print(f"-b selected, loading {intree} with quoted node values.")
#     try:  # check if intree is Newick
#         Tree(intree, format=1)  # check if this works
#     except FileNotFoundError:  # is this try/except loop necessary?
#         raise
# else:
#     try:  # check if intree is Newick
#         Tree(intree)
#     except FileNotFoundError:
#         raise

try:
    t = Tree(intree)
except FileNotFoundError:
    sys.exit(f"File {intree} not found. Exiting.")
except NewickError:
    print(f"{intree} contains unusual boostraps, loading as strings (format=1):")
    try:
        t = Tree(intree, format=1)
        weirdboots = True
    except NewickError:
        sys.exit("Unable to parse tree file with quoted strings(format=1) either, is it in Newick?")
        # if loaded as list, skip offending tree(s)
else:
    weirdboots = False


# Log loaded files:
# including args settings, eg. overwrite!

def get_time():
    return time.strftime("%d %b %Y, %X")


write_to_logfile(["", "", f"Generated by {sys.argv[0]} run on {get_time()}.",
                  f"Input tree is: {intree}", f"File for full taxon names: {tax_to_fullname}",
                  f"Taxonomic groups specified in: {taxgroups}", f"Colours specified in: {colourlist}",
                  f"Weird boots: {weirdboots}, Overwrite: {args.overwrite}, Coverage stats: {args.cov_stats}",
                  f"Outfile: {args.outfile}, Log file: {logfile}", ""])


#####################
# Parse input files #
#####################

#get unique values from nested dictionary
def get_unique_values(dic, key2):
    s=[]
    for key1 in dic:
        s.append(dic[key1][key2])
    return(set(s))


# Load colour palette
Colours = {}
with open(colourlist, "r") as f:
    for l in csv.reader(f, delimiter='\t'):
        Colours[l[0]] = l[1]


# load names

Fullnames = {}
with open(tax_to_fullname, "r") as f:
    for l in csv.reader(f, delimiter='\t'):
        Fullnames[l[0]] = l[1]  # change to l[2] for complete species names if needed

# print(Fullnames)  # add this to log file

# MakeDictionary
TaxInfo = {}
with open(taxgroups, "r") as f:
    for row in csv.reader(f, delimiter='\t'):
        TaxInfo[row[0]] = {"BaseName": row[0],
                           "FullName": Fullnames[row[0]],
                           "Group": row[1],
                           "Colour": Colours[row[1]]
        }
#       TaxNames[l[0]] = l[1:] # to read the rest


# figure out a more sensible way to log stuff
# with open(logfile, "a") as log:  # refactor with log= instead
#     log.write("Parsed information is:\n")
#     for i in TaxInfo:
#         log.write(i+":\t"+TaxInfo[i]["FullName"]+"\t"+TaxInfo[i]["Group"]+"\t"+TaxInfo[i]["Colour"]+"\n")

TaxInfo_log = ["Parsed information is:"]
for i in TaxInfo:
    TaxInfo_log.append(i+":\t"+TaxInfo[i]["FullName"]+"\t"+TaxInfo[i]["Group"]+"\t"+TaxInfo[i]["Colour"])
write_to_logfile(TaxInfo_log)

# load coverage stats
def parse_cov_stats(file, genes_sites):  # could also run check_cov_stats inside this function
    """
    Parses gene and sites coverage stats file.
    :param file: tab-delimited coverage stats file as above
    :param genes_sites: output from check_cov_stats as list of [True, True], etc
    :return: Dictionary formatted as TaxonName : gene or perc sites
    """
    cov_dic = {}
    with open(file, "r") as cov_f:
        if genes_sites == [True, True]:  # read both genes and sites
            for l in csv.reader(cov_f, delimiter='\t'):
                try:
                    cov_dic[l[0]] = {"GeneCount": int(l[1]), "PercSites": float(l[2])}
                except ValueError:
                    print(f"Non-numeric value (or non-whole-number value for gene count) detected, skipping coverage stats!\
                          \noffending line: {l}")
                    return None
            return cov_dic
        elif genes_sites == [True, False]:  # genes only
            for l in csv.reader(cov_f, delimiter='\t'):
                try:
                    cov_dic[l[0]] = {"GeneCount": int(l[1])}
                except ValueError:
                    print(f"Non-whole-number value for gene count detected, skipping coverage stats! \noffending line: {l}")
                    return None
            return cov_dic
        elif genes_sites == [False, True]:  # sites only
            for l in csv.reader(cov_f, delimiter='\t'):
                try:
                    cov_dic[l[0]] = {"PercSites": float(l[2])}
                except ValueError:
                    print(f"Non-numeric value for percent sites detected, skipping coverage stats! \noffending line: {l}")
                    return None
            return cov_dic
        else:
            print(f"Something is wrong, skipping coverage stats!")
            return None


def test_cov_values(cov_stats_dic):  # appropriate to use void function?
    """
    Tests loaded values in coverage stats for gene counts over specified maximum and non-percentages (note: will not catch
percentages written as decimals, would appear as extremely low values instead.
    :param cov_stats_dic: coverage stats dictionary generated by parse_cov_stats
    :return: void
    """
    try:
        maxgene = max(get_unique_values(cov_stats_dic, "GeneCount"))
        if maxgene > cov_stats_maxgenes:
            print(f"Found gene count({maxgene}) exceeding specified max genes({cov_stats_maxgenes}), skipping coverage stats.")
            args.cov_stats = None
    except KeyError:  # could be intentionally missing
        pass
    try:
        maxperc = max(get_unique_values(cov_stats_dic, "PercSites"))
        if maxperc > 100:
            print(f"Percentages must be between 0 and 100; value >100 detected ({maxperc}), skipping coverage stats.")
            args.cov_stats = None
    except KeyError:  # could be intentionally missing
        pass  # this might also be a dumb way to do this


if args.cov_stats:
    CovStats = parse_cov_stats(cov_stats_file, incl_genes_sites)
    # test coverage stats
    if CovStats == None:
        args.cov_stats = None
    else:
        test_cov_values(CovStats)  # to check blatant errors in coverage stats file


#################################################
# Load tree and determine taxon colouring order #
#################################################

# Load tree

# t = Tree(intree) replaced above!
leaflist = t.get_leaf_names()  # list of taxa in present tree


# Set outgroup ## Make this user-defined later
def set_outgroup(tree, tax1, tax2):
    ancestor = tree.get_common_ancestor(tax1, tax2)
    tree.set_outgroup(ancestor)


try:
    set_outgroup(t, 'Bodosalt', 'TriPCT')
except ValueError:
    print("\nOne of outgroup defining taxon pairs not found in tree, using midpoint rooting!\n")
    t.set_outgroup(t.get_midpoint_outgroup())


# Group [Taxon: group] by taxon (later change input file/how it's parsed)
def group_by_taxon_group(dic):
    d = {}
    for k in dic:
        if dic[k] not in d:
            d[dic[k]] = [k]
        else:
            d[dic[k]].append(k)
    return d
# format of 'rearranged' dictionary: {organism : group}


# first, extract all group info from TaxInfo:  ## probably don't need
# def get_unique_values(dic, key2):
#     s=[]
#     for key1 in dic:
#         s.append(dic[key1][key2])
#     return(set(s))
# for i in get_unique_values(TaxInfo, "Group"):


# extract from TaxInfo a {taxon, groupname}
rearranged = {}
for key in TaxInfo:
    rearranged[key] = TaxInfo[key]['Group']
TaxaListedbyGroups = group_by_taxon_group(rearranged)

logtext = []  # list of strings to print to logfile
# check against taxa in current tree, remove those that are missing, and remove blank groups (print them first)
logtext.append(f"Trimming source list for taxa not in tree: {intree}")
AbsentGroups = []
for groupname in TaxaListedbyGroups:
    # remove taxa that are not in the tree being processed; if later versions use multiple trees
    # with one template, make sure the TaxaListedbyGroups value is regenerated for each tree!
    # Must happen before rulepairs generation -- more elegant way would be to have rule generation
    # account for missing taxa and groups.
    PresentTaxa = [name for name in TaxaListedbyGroups[groupname] if name in leaflist]
    if len(PresentTaxa) == 0:
        logtext.append(f"{groupname} entirely absent from {intree}.")
        print(logtext[len(logtext) - 1])
        AbsentGroups.append(groupname)
        pass
    elif len(PresentTaxa) != len(TaxaListedbyGroups[groupname]):
        logtext.append("Missing taxa for %s: %s" %
                       (groupname, " ".join([n for n in TaxaListedbyGroups[groupname] if n not in PresentTaxa])))
        # print(logtext[len(logtext) - 1])
        TaxaListedbyGroups[groupname] = PresentTaxa
for groupname in AbsentGroups:
    TaxaListedbyGroups[groupname].pop()

write_to_logfile(logtext)
logtext = []

# print(TaxaListedbyGroups)  # log this?


# ## Order user-defined groups so that paraphyletic taxa don't colour branches after nested clades are done
# assumes monophyly of groups (incl paraphyly, but NOT polyphyly!)
# get corresponding tree
def get_subtrees(grouplist):
    # expected input: {groupname : [taxon1, taxon2, taxon3...]}
    treedic = {}
    for k in grouplist:
        treedic[k] = t.get_common_ancestor(grouplist[k])
    return treedic


def order_groups(dic):
    """
    Orders taxon groups in nesting order to colour innermost groups last in given tree
    :param dic: expected input: {groupname : subtree}
    :return: ordered list of rule pairs, nesting group containing nested in that order
    """
    ordered_list = []
    for groupname in dic.keys():
        for othergroupname in dic:
            if groupname != othergroupname:
                if dic[groupname] in dic[othergroupname].traverse():
                    logtext.append(f"{groupname} is contained by {othergroupname}")
                    ordered_list.append([othergroupname, groupname])  # in that order, nesting group contains nested
    return ordered_list


rulepairs = order_groups(get_subtrees(TaxaListedbyGroups))
# later: add detection for polyphyletic groups, to colour separately, somehow
write_to_logfile(logtext)


# to build longest paths of rule pairs
def flatten_rules(rulesets):
    #  needs several iterations if multiple depths
    for i, v in enumerate(rulesets):
        for o in rulesets:
            if v[(len(v)-1)] == o[0]:
                rulesets[i] = list(itertools.chain(*[v, o[1:]]))


# to find longest rule (this returns length as int)
def max_length(list0):
    longest = max(len(elem) for elem in list0)
    return longest


# to then find rule sets by length:
def get_element_by_length(list0, length):
    outlist = []
    for i in list0:
        if len(i) == length:
            outlist.append(i)
    return outlist


# then clean up short bits (this is like genome assembly, sort of)
def simplify_rules(longelement, rulesets):
    for v in rulesets:
        if v != longelement:
            if all(elem in longelement for elem in v):
                rulesets.remove(v)


# this cleans up the rules non-repeating elements only:
l = max_length(rulepairs)
ordered_rule_list = []
while l > 1:
    for x in get_element_by_length(rulepairs, l):
        simplify_rules(x, rulepairs)
        for groupname in x:
            if groupname not in ordered_rule_list:
                ordered_rule_list.append(groupname)
        # to then add to ordered_list -- this is an untested approach!
    l -= 1
# now the rest of the groupnames:


for groupname in TaxaListedbyGroups.keys():
    if groupname not in ordered_rule_list:
        ordered_rule_list.append(groupname)

# then refer to ordered_rule_list in colouring step!


#############################
# Colouring/formatting tree #
#############################

if weirdboots:
    # first value is assumed to be bootstrap!
    for node in t.traverse():
        if not node.is_leaf():
            if node.name != "":
                if len(node.name.split("/")) > 2:
                    print("More than two node values detected, treating all as name (may add support for more values later)")
                elif len(node.name.split("/")) == 1:
                    print("Not separable by /, treating all as one name.")

                # insert problematic tree filename (exactly as called in script) here to use second value after /
                # elif intree == "Trees/66t254gh05g02FINAL.EsModel.fasta.treefile":  # very specific problem, delete after use
                #     node.support = float(node.name.split("/")[1])
                #     node.name = ""

                else:
                    node.support = float(node.name.split("/")[0])
                    # node.name = node.name.split("/")[1]  # normal mode, uncomment this
                    node.name = round(float(node.name.split("/")[1]), 1)  # for gCF, comment away after

            else:
                pass


# Add full names to tree:
for leaf in t:
    leaf.add_features(FullName=str("  "+Fullnames[leaf.name]))  # make this optional if not args.tax_to_fullname
    if args.cov_stats is not None:  # make this work without input for cov_stats
        try:
            leaf.add_features(CovGenes=CovStats[leaf.name]["GeneCount"])
        except (KeyError, NameError) as er:
            # print(er)
            pass
        try:
            leaf.add_features(CovSites=CovStats[leaf.name]["PercSites"])
        except (KeyError, NameError) as er:
            # print(er)
            pass


def colour_node(node, group):
    if node in t.get_common_ancestor(TaxaListedbyGroups[group]).traverse():
        node.img_style["vt_line_color"] = Colours[group]
        node.img_style["hz_line_color"] = Colours[group]
        node.img_style["fgcolor"] = Colours[group]
        # add default value: NestedDic[key1].get(key2, DefaultValue)

######################################
# Move this somewhere else (or delete) -- this is playing around with gCF
# def gCF_min_max(t):
#     boots = []  # dumb var name, rename
#     for n in t.traverse():
#         if not n.is_leaf() and n.name != "":
#             boots.append(float(n.name))
#     min_boots = min(boots)
#     max_boots = max(boots)
#     return [min_boots, max_boots]
#
# import numpy  # remove this later
# boots_min_max = gCF_min_max(t)
# boots_min_max = numpy.cbrt(gCF_min_max(t))
######################################

# Colouring, removing nodes, etc
def my_layout(node):
    node.img_style["hz_line_color"] = '#797b7d'  # change to refer to DefaultSettings
    node.img_style["vt_line_color"] = '#797b7d'
    node.img_style["hz_line_width"] = 2  # old setting: 1
    node.img_style["vt_line_width"] = 2
    node.img_style["size"] = 0
    for group in ordered_rule_list:
        colour_node(node, group)
    if not node.is_leaf() and weirdboots:  # for GCF and other additional support values
        extra_node_label = TextFace(str(node.name), fsize=(tree_leaf_font_size),
                                 fstyle="bold",
                                 ftype="Myriad pro",
                                 fgcolor="#0000D5"
                                 )
        extra_node_label.margin_bottom = 8  # def 12; note, if want to appear below bootstraps,
        # change order with bootstrap block
        parent_branch_midpoint = (node.dist / 2) * ts.scale  # later use this for .jplace placements!
        # print(f"branch midpoint {parent_branch_midpoint}")
        if parent_branch_midpoint > 12:  # moves text to middle of long branches; looks weird for short branches
            extra_node_label.margin_right = parent_branch_midpoint
        node.add_face(extra_node_label, column=0, position="float")

        # # # GCF # # #
        #
        # testing branch-thickness idea!
        # if node.name == "":
        #     pass
        # else:
        #     max_branch_thickness = 8
        #     node_gCF = float(node.name)
        #     ratio = numpy.cbrt(node_gCF / boots_min_max[1])
        #     # gCF_width = ratio * max_branch_thickness
        #     # gCF_width = numpy.sqrt(ratio * max_branch_thickness)
        #     gCF_width = ratio * max_branch_thickness
        #     node.img_style["hz_line_width"] = gCF_width
        #     node.img_style["vt_line_width"] = 1
        #
        #
        #

    if not node.is_leaf() and node.support > 99.9:
        node.img_style["size"] = 5  # previously: 4
        # node.img_style["size"] = gCF_width + 3  # for messing around with GCF branch thickness stuff (below)
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = node.img_style["hz_line_color"]  # for default nodes with full support
    else:
        if not node.is_leaf():
            node.img_style["size"] = 0
            if node.support > 50:  # set to 0 if show all branches
                node_text = TextFace(str(int(node.support)), fsize=(tree_leaf_font_size - 2),
                                     fstyle="bold",  # bold is not an option
                                     ftype="Myriad pro",  # Verdana, Arial, Myriad
                                     # fgcolor=node.img_style["hz_line_color"]  # to match text colour to branch colour
                                     )
                # node_text.margin_bottom = 12  # make function of leaf font size, branch thickness (2+1+font size?)
                node_text.margin_top = 6  # def: 12 switch to this to display numbers below branches
                # node_text.opacity = 0.5
                node.add_face(node_text,
                              column=0,
                              position="float",  # avoids branch length extension issue (with branch-top)!
                              )
            else:
                pass
        else:
            node.img_style["size"] = 0
            leaf_face = AttrFace(attr="FullName", fsize=tree_leaf_font_size, fstyle="italic",
                                 # text_prefix=" "  # remove text_prefix!!!
                                 )
            faces.add_face_to_node(face=leaf_face, node=node, column=1, position="branch-right")
            if args.cov_stats:
                try:
                    cov_bar1 = render_bar_graph(value=node.CovGenes,
                                                max_value=cov_stats_maxgenes,
                                                colour=node.img_style["hz_line_color"],
                                                border_colour=None,  # node.img_style["hz_line_color"]
                                                bar_height=3,
                                                max_width=100)

                # cov_bar2 = render_bar_graph(value=node.CovSites,  # to show max value  # to represent by single graph
                #                             max_value=100,
                #                             colour=node.img_style["hz_line_color"],
                #                             border_colour="#FFFFFF",
                #                             max_width=100)

                    cov_bar2 = render_perc_bar_graph(node.CovSites,
                                                     fore_colour=node.img_style["hz_line_color"],
                                                     bar_height=3)

                    cov_bar1.margin_left = 35  # edit this to auto-detect second/third longest branches + longest text
                    cov_bar2.margin_left = 35  # edit this to auto-detect second/third longest branches + longest text
                    cov_bar2.opacity = 0.5  # make second bar lighter
                    cov_bar1.margin_bottom = 0.5
                    cov_bar2.margin_top = 0.5
                    faces.add_face_to_node(cov_bar1, node=node, column=3, position="branch-right", aligned=True)  # can add border here
                    faces.add_face_to_node(cov_bar2, node=node, column=3, position="branch-right", aligned=True)

                except (KeyError, NameError) as er:
                    write_to_logfile(er)

            # GLITCH: uses longest branch + taxon for scaling; somehow tie in scaling factor to add padding for bar graph
            # print(node.faces)  # for debugging


def render_bar_graph(value, max_value, colour, max_width=100, bar_height=5, border_colour="#FFFFFF"):
    """
    Generates bar object for adding coverage stats to the right of the tree
    :param value: value to depict in this instance (eg. node.CovGenes or node.CovSites)
    :param max_value: maximum value for feature; 100 for percentages, cov_stats_maxgenes for gene count.
    :param colour: either node.img_style['hz_line_color'] or hex value of choice, eg #000000 for black
    :param max_width: max width of bars in pixels
    :param bar_height: thickness of bars, in pixels. 5 (default) for single, use 2 for double
    :param border_colour: colour of bar graph border; default = #FFFFFF; set to bar colour when comparing w stacked bar
    :return: RectFace object for particular node
    """
    width_factor = max_width / max_value
    cov_bar = RectFace(height=bar_height,
                       width=(value * width_factor),  # height = 5 IF only one feature!
                       fgcolor=border_colour,  # white border, or user-specified
                       bgcolor=colour)
    return cov_bar


def render_perc_bar_graph(value, max_value=100, fore_colour="#B8B8B8", back_colour="#D4D4D4",
                          max_width=100, bar_height=5, border_colour="#FFFFFF"):
    """
    Generates stacked bar object for adding coverage stats to the right of the tree
    :param value: input value (eg. node.CovGenes or node.CovSites)
    :param fore_colour: grey; can also be node colour
    :param max_value: default 100 for percentages.
    :param back_colour: background colour, default light grey (for max - value)
    :param max_width: width of total bar object, in pixels
    :param bar_height: height (thickness) of total bar object, in pixels
    :param border_colour: colour of border around total object, default white.
    :return: StackedBarFace object for particular node
    """
    cov_bar = StackedBarFace([value, (max_value-value)],
                             colors=[fore_colour, back_colour],  #F0F0F0 too light
                             height=bar_height,
                             width=max_width,
                             line_color=border_colour)  # white border, or user-specified
    return cov_bar
# test = StackedBarFace([55,45],width=100,height=5, colors=["#F0F0F0","#D4D4D4"], line_color="#FFFFFF")


def find_scaling_factor(tree, paper_format="letter", leaf_font_size=10):
    """
    :param tree: input tree
    :param paper_format: 'letter', 'A4', custom width/height ratio, 'letter-landscape', 'A4-landscape'
    :param leaf_font_size: font size of leaf text, in pt (default = 10)
    :return: scaling factor for ts.scale
    """
    numtax = len(tree.get_leaf_names())
    max_branch = tree.get_farthest_leaf()[1]  # second element is branch length to farthest leaf
    if paper_format == "letter":
        paper_width_to_height = 0.772
    elif paper_format == "A4":
        paper_width_to_height = 0.707
    elif paper_format == "letter-landscape":
        paper_width_to_height = 1.294
    elif paper_format == "A4-landscape":
        paper_width_to_height = 1.414
    else:
        paper_width_to_height = int(paper_format)
    scale = numtax * (4 / 3) * (leaf_font_size + 1) * paper_width_to_height / max_branch
    return scale


def set_scalebar_length(tree, ratio_max_length=0.20):
    """
    Determines optimal scalebar length for ts.scale_length
    :param tree: input tree
    :param ratio_max_length: percent of longest branch to set scalebar length to, 20% default. 0.1 for 10%, etc.
    :return: Optimised scale setting
    """
    max_branch = tree.get_farthest_leaf()[1]
    if max_branch * ratio_max_length >= 1:
        rounding_digits = 0
    elif max_branch * ratio_max_length >= 0.1:
        rounding_digits = 1
    elif max_branch * ratio_max_length >= 0.01:
        rounding_digits = 2
    elif max_branch * ratio_max_length >= 0.001:
        rounding_digits = 3
    else:
        rounding_digits = 4
    return round(max_branch * ratio_max_length, rounding_digits)
    # update for next time: convert output value to string and count 0s after . until first non-0 --> round(ndigits)


if args.cov_stats:
    print(f"-s selected, making bar graph of gene and site coverage \
specified in {cov_stats_file}, total of {cov_stats_maxgenes} genes expected")


ts = TreeStyle()
ts.layout_fn = my_layout
ts.optimal_scale_level = "full"  # other option:  med, full  FULL is better!
ts.show_branch_support = False
ts.show_leaf_name = False

ts.scale = find_scaling_factor(t, paper_format="letter", leaf_font_size=10)

ts.show_scale = True
ts.scale_length = set_scalebar_length(t, ratio_max_length=0.10)
print(f"scale is {ts.scale}")

# add later: optimise spacing via ts.branch_vertical_margin and font size; choose paper size,
# auto scale (full or mid) for exceptionally weird topologies/too many taxa/etc

# ts.branch_vertical_margin = 8    # 8 with 66 taxa at 10pt font with 1.0090619846 maxbranch
# ts.scale = numtax * (4/3) * (tree_leaf_font_size+1) * letter_width_to_height / max_branch

"""
numtax *   # number of taxa in tree
(4/3) *   # scaling factor for pixels/point
(tree_leaf_font_size+1) *   # font size + 1 (default spacing margin, points)
 letter_width_to_height  # page format w/h
 / max_branch  # maximum branch length
"""

# ts.legend.add_face(TextFace("0.5 support"), column=1)
# t.show(tree_style=ts)
# If tree shown AND rendered, bootstrap support values appear twice on top of each other. Do one OR other.


################
# Render tree! #
################

# for incremental file naming of rendered output:
def extract_number_from_outfile(f):
    s = re.findall("\d+.svg", f)
    return int(s[0].split(".")[0]) if s else -1, f  # returns -1 if no number found at end of file


def generate_next_filename(infile, extension):
    """
    Generates new filename incremented by 1; creates *1.{extension} if extant outfile unnumbered
    :param infile: input tree file, as called in script; path from current working directory will be prepended.
    :param extension: file extension; will be simply added to the end of infile+number
    :return: new filename as {path/to/current/working/directory}/{infile}{1,2,3,10,...etc}.{extension}
    """
    cwd = getcwd()
    existing_svgfiles = filter(path.isfile, glob.glob(f"{cwd}/{infile}*.{extension}"))
    fn = max(existing_svgfiles, key=extract_number_from_outfile)  # find highest number
    fdir = path.split(fn)[0] + "/"  # extract directory/path to filename
    fn = path.split(fn)[1]  # extract filename
    try:
        fn_int = re.findall(f"\d+.{extension}", fn)[0]  # get filename before .svg (or extension)
    except IndexError:
        print("Starting new numbered file")  # kludge for when file exists but not yet numbered
        # fn_int = "0.svg"  # this will become 1.svg
        # kludge for when file exists but not yet numbered
        fn_int = f"0.{extension}"  # this will become 1.svg
        fn = ".".join(fn.split(".")[0:(len(fn.split("."))-1)])+fn_int  # -1 to trim off svg
        # there is certainly a less idiotic way to do this
    fn_int_next = str(int(fn_int.split(".")[0]) + 1) + f".{extension}"  # increment number in file
    new_name = fdir + re.sub(str(fn_int), fn_int_next, fn)  # substitution
    return new_name


# intree with path  # this can be written better as one/two functions
if args.outfile:  # doing it this way will break when called from another script; refactor later
    if args.overwrite or not path.isfile(args.outfile):
        outfile_name = args.outfile
    else:
        exit(f"Error: {args.outfile} already exists! Quitting.")
        # error+quit instead of saving under default name prob better for pipelines
else:
    if not args.overwrite:
        if path.isfile(f"{intree}.svg"):
            outfile_name = generate_next_filename(intree, "svg")
        else:
            outfile_name = f"{intree}.svg"
    else:
        print(f"Overwrite mode on.")
        outfile_name = f"{intree}.svg"

t.render(outfile_name,
         tree_style=ts,
         w=800,
         #  w=ts.scale * 0.8,  # seems to render best at 400 px or mm; 375-425 tolerable  #set to width in scale?
         units="px")

write_to_logfile([f"Rendered file saved to: {outfile_name}"])  # move to logfile
print(f"Rendered file saved as: {outfile_name}")

##################################
# Move to associated notes file! #
##################################

# # Save only one render per run; otherwise, errors accumulate (weird defect around support val when rendered twice)
# t.render(f"{intree}1.svg",  #set up for overwrite vs. increment file number
#          tree_style=ts,
#          #w=300,  # do not make much bigger than 400 for svg rendering, the renderer starts to suck
#          # 150 -- much too small!
#          # a difference between rendering of SVG in Preview vs. Illustrator...!?
#          w=400,
#          units="mm"
#          # dpi=
#          )

# pdf text rendering defective, each letter as independent object

# t.render(f"{intree}.pdf",
#          # optimal_scale_level="full"
#          tree_style=ts,
#          w=800,
#          units="mm"
#          # dpi=
#          )
