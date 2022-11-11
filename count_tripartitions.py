#!/usr/bin/env python

from ete3 import Tree
import argparse

# make command line interface
parser = argparse.ArgumentParser(
    description="""Count tripartitions in a tree or set of bootstrap trees."""
)
parser.add_argument(
    "-t", 
    dest="tree_file",
    metavar="TREE_FILE",
    # required=True, 
    help="Tree or trees in Newick format"
)
args = parser.parse_args()


# define count dict
count = {}


# define tree_total count
total_tree_count = 0


# open tree file
with open(args.tree_file, "r") as f:

	# each line is a bootstrap / mcmc tree
	for line in f:

        # update total tree count
		total_tree_count += 1

    	# read in tree
		t = Tree(line)
                
		# cache tree
		## keys = node objects, values = sets of leaf names
		cache = t.get_cached_content(store_attr="name")
		# store all  leave names in set
		all_leaves = set( [leaf.name for leaf in t.iter_leaves()] )

		# traverse nodes in bootstrap tree
		for node in t.traverse("postorder"):
			
			# skip if leaf node
			if node.is_leaf():
				continue

			# define 3 sets of taxa that emerge from node
			## taxa associated with branch one:
			children_one = cache.get(node.children[0])
			## taxa associated with branch two:
			children_two = cache.get(node.children[1])
			## all other taxa
			all_other_children = all_leaves.difference(children_one.union(children_two))

			# define tripartition as three strings of taxa
			tripartition = "____".join(sorted( 
				[ 
					",".join(sorted(children_one)), 
					",".join(sorted(children_two)), 
					",".join(sorted(all_other_children)) 
				] 
			))
			# print(type(tripartition))
			# print(tripartition)

			# add one to tripartition count
			if tripartition in count:
				count[tripartition] += 1
			else:
				count[tripartition] = 1


# print output
for tripartition,count in count.items():
	percent_trees = round(count / total_tree_count * 100)
	print(str(percent_trees) + "\t" + tripartition)
