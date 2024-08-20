#!/usr/bin/env python

"""
Copyright 2023 Joran Martijn.

The development of this and other bioinformatic tools in the Roger lab was funded by
Discovery grant RGPIN-2022-05430 from the Natural Sciences and Engineering Research Council of Canada
awarded to Andrew J. Roger.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.
If not, see <https://www.gnu.org/licenses/>.
"""

# default libraries
import argparse

# stuff you need to install
import gffutils
import re

# NOTE:
#   The script assumes that exon features have an ID 
#   that includes the word 'exon'

# make command line interface
parser = argparse.ArgumentParser(description="Add intron features to a GFF3 file")

parser.add_argument(
    "-g", 
    dest="gff3_file",
    metavar="GFF3_FILE",
    required=True, 
    help="Input GFF3 file"
)

args = parser.parse_args()


def main(args):

    # load gff file into memory
    db = gffutils.create_db(
        data=args.gff3_file,
        dbfn=':memory:',
        merge_strategy="create_unique"
    )

    # infer and add intron features
    introns = list(db.create_introns())

    # edit ID attributes of introns
    ## created intron features have ID=exon1,exon2 as format by default
    ## 'exon1' and 'exon2' come from the attribute IDs from the exon features
    ## this breaks .update(), so we need to adjust their IDs

    ## the order of listed introns is alphabetical sort, so
    ## exon7,exon8
    ## exon8,exon9
    ## exon10,exon9 !!
    ## exon10,exon11

    prev_mrna = ''

    for i in introns:

        # if exons don't have IDs, resultant introns won't either,
        # so we'll have to create them here
        if not 'ID' in i.attributes:
            mrna = i.attributes['Parent'][0]
            if prev_mrna == mrna:
                n += 1
            else:
                prev_mrna = mrna
                n = 1
            i.attributes['ID'] = [f'{mrna}.exon{n}',f'{mrna}.exon{n+1}']

        # print(i)

        # the edits below only work with gffutils 0.11!
        # in gffutils 0.11, i.attributes['ID'] = ['exon1','exon2']
        # in gffutils 0.12, i.attributes['ID'] = ['exon1-exon2']
        # so, if running with 0.12, comment out the next 3 code lines
        # sort exon10,exon9 -> exon9,exon10
        sorted_exon_names = sorted( i.attributes['ID'], key=natural_sort_key )

        # remove the second exon of the ID
        i.attributes['ID'].pop(1)

        # give the intron feature an intron ID!
        i.attributes['ID'][0] = sorted_exon_names[0].replace('exon','intron')

        # change source field from 'gffutils_derived' to 'gffutils'
        i.source = 'gffutils'

    # update db with new introns
    db.update(introns)

    # print updated gff record in logical order
    for g in db.features_of_type('gene', order_by=('seqid', 'start')):
        print()
        print(g)
        for m in db.children(g, featuretype=('mRNA','tRNA')):
            print(m)
            for f in db.children(m, order_by='start'):
                print(f)


# we need a natural sort function
# to sort exon10,exon9 -> exon9,exon10
# I don't understand it, got it from StackOverflow
def natural_sort_key(s):
    _nsre = re.compile('([0-9]+)')
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]


if __name__ == '__main__':
    main(args)
