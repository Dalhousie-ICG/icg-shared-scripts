#!/usr/bin/env python

import argparse
import csv


# make command line interface
parser = argparse.ArgumentParser(
    description="""
				Collect tripartition counts of tripartitions
				of interest in multiple tripartition files
				"""
)
parser.add_argument(
    "-t", 
    dest="tripartitions",
    metavar="TRIPARTITIONS",
    type=str,
    nargs="*", # setting nargs to * will allow -t to take any number of arguments
    required=True, 
    help="Counts of all tripartitions"
)
parser.add_argument(
	"-m",
	dest="map",
	metavar="TRIPARTITIONS_OF_INTEREST",
	required=True,
	help="Tripartions we wish to track"
)
args = parser.parse_args()


# store tripartitions of interest
with open(args.map, "r") as f:
	map = {}
	for label, tripartition in csv.reader(f, delimiter="\t"):
		# ignore commented labels
		if not label.startswith('#'):
			map[tripartition] = label


# store list of files
files = args.tripartitions


# set counts to 0 per tripartition file per label
## so 
## {
##  'file1': {'label1': 0, 'label2': 0, etc},
##  'file2': {'label1': 0, 'label2': 0, etc},
##  'etc': {etc} 
## }
tracked_counts = {}
for file in files:
	tracked_counts[file] = {}
	for label in map.values():
		tracked_counts[file][label] = 0
# print(tracked_counts)


# parse files and read tripartition counts into tracked_counts dict
for file in files:
	with open(file, "r") as f:
		for count, tripartition in csv.reader(f, delimiter="\t"):
			# check if tripartition in this line is one we want
			if tripartition in map:
				# retrieve label (e.g. 'alpha_sister') from map dict
				tracked_label = map.get(tripartition)
				# link label to tripartition count
				## (overwrite 0 if label was found)
				tracked_counts[file][tracked_label] = int(count)


# print counts of interest in human readable table
labels = map.values()
print( "file", "\t", "\t".join(labels) )
for file in files:
	print(file, end='')
	for label in labels:
		print("\t", tracked_counts[file][label], end='')
	print(end='\n')
