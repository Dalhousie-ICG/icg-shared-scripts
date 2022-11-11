#!/usr/bin/env python

import argparse
import gffutils

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

# load gff file into memory
db = gffutils.create_db(
    data=args.gff3_file,
    dbfn=':memory:',
    merge_strategy="error"
)

# infer and add intron features
introns = list(db.create_introns())

# edit ID attributes of introns
## created intron features have ID=exon1,exon2 as format by default
## this breakse .update(), so we need to adjust their IDs
for i in introns:

    # remove the second exon of the ID
    i.attributes['ID'].pop(1)

    # give the intron feature an intron ID!
    i.attributes['ID'][0] = i.attributes['ID'][0].replace('exon','intron')

    # change source field from 'gffutils_derived'
    # to 'gffutils'
    i.source = 'gffutils'

# update db with new introns
db.update(introns)

# print updated gff record in logical order
for f in db.all_features(order_by=('seqid','start','attributes')):#, limit='tig00000012:1-400000'):
    if f.featuretype == 'gene':
        print()
        print(f)
    else:
        print(f)
