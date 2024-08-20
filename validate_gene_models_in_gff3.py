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

# standard modules
import argparse
import sys

# stuff you need to install
import gffutils
import pyfaidx
from Bio.Seq import Seq
import regex as re


# argument parser
parser = argparse.ArgumentParser(
            description=
        '''
        Do some common sense sanity checks on the gene models.

        Does the gene feature have any children features?
        Do gene coordinates match with mRNA & tRNA coordinates?
        Does the start coordinate of the first CDS/exon match with the gene start?
        Does the end coordinate of the last CDS/exon match with the gene end?
        Is the cumulative length of all exons a multiple of 3? 
        Does the translate CDS start with an M and end with a *?
        Are there any premature * in the CDS?

        Optionally print the CDS sequences in FASTA format.
        Genes with multiple CDSs/exons are merged and then translated into a single protein sequence.
        ''',
        formatter_class=argparse.RawTextHelpFormatter
        )
parser.add_argument(
        "-g", "--gff3",
        dest='gff3_file',
        type=str,
        required=True,
        help="Input GFF3 genome file")
parser.add_argument(
        "-f", "--fasta",
        dest='fasta_file',
        type=str,
        required=True,
        help="Input FASTA file")
parser.add_argument(
        "-b", "--blasto",
        dest='is_blastocystis',
        action='store_true',
        required=False,
        help="Toggles checking for T....TGTTTGTT motif")
parser.add_argument(
        "--print_proteins",
        action='store_true',
        help="Toggles printing of protein CDS sequences in FASTA format")
args = parser.parse_args()


# load GFF3 and FASTA files into memory
db = gffutils.create_db(args.gff3_file, dbfn=':memory:', merge_strategy='create_unique')
fasta = pyfaidx.Fasta(args.fasta_file)

errors = []

# one CDS sequence per gene
for gene in db.features_of_type('gene', order_by=('seqid', 'start')):
    # print()
    # print(gene.id)

    # check whether the gene has any children
    # e.g., 'tRNA' or 'mRNA' or 'CDS' features
    try:
        child = next( db.children(gene) )
    except StopIteration:
        print(f"Gene {gene.id} on {gene.seqid} lacks children features!")

    # store all featuretypes associated with gene in memory
    gene_features = [ child.featuretype for child in db.children(gene) ]
    # print(gene_features)


    # check if there are any RNA children gene features
    if all(rna not in gene_features for rna in ['mRNA', 'tRNA', 'rRNA']):
        warning_message = (
            f"Warning: Gene {gene.id} on {gene.seqid} lacks RNA features!\n"
            "This in principle OK but doesn't represent biology very well"
        )
        # print('', file=sys.stderr)
        # print(warning_message, file=sys.stderr)
        # print('', file=sys.stderr)

    # if no RNAs, are 'CDS' the only featuretype among children?
    if set(gene_features) == {'CDS'}:
        warning_message = (
            f"Warning: Gene {gene.id} on {gene.seqid} has only CDS features!\n"
            "This in principle OK but doesn't represent biology very well"
        )
        # print(warning_message, file=sys.stderr)


    # checks for tRNA genes
    if 'tRNA' in gene_features:

        # check if gene and tRNA coordinates match
        trna = next( db.children(gene, featuretype='tRNA') )
        if trna.start != gene.start:
            errors.append( ValueError(f"Gene {gene.id} on {gene.seqid} has tRNA feature but start coordinates do not match") )
        elif trna.end != gene.end:
            errors.append( ValueError(f"Gene {gene.id} on {gene.seqid} has tRNA feature but end coordinates do not match") )

        # check exon feature coordinates
        exon = next( db.children(gene, featuretype='exon') )
        if exon.start != gene.start:
            errors.append( ValueError(f"Gene {gene.id} on {gene.seqid} has exon feature but start coordinates do not match") )
        elif exon.end != gene.end:
            errors.append( ValueError(f"Gene {gene.id} on {gene.seqid} has exon feature but end coordinates do not match") )

        # move to next gene
        continue


    # checks for mRNA genes
    elif 'mRNA' in gene_features:

        # do gene and mRNA features coordinates match?
        mrna = next( db.children(gene, featuretype='mRNA') )
        if mrna.start != gene.start:
            errors.append( ValueError(f"Gene {gene.id} on {gene.seqid} has mRNA feature but start coordinates do not match") )
        elif mrna.end != gene.end:
            errors.append( ValueError(f"Gene {gene.id} on {gene.seqid} has mRNA feature but end coordinates do not match") )

        # is there an equal number of CDS and exon features?
        exon_count, CDS_count = gene_features.count('exon'), gene_features.count('CDS')
        if exon_count > CDS_count:
            errors.append( ValueError(f"Gene {gene.id} on contig {gene.seqid} has more exon features ({exon_count}) than CDS features ({CDS_count})") )
            continue
        elif exon_count < CDS_count:
            errors.append( ValueError(f"Gene {gene.id} on contig {gene.seqid} has more CDS features ({CDS_count}) than exon features ({exon_count})") )
            continue

        ## exon boundaries should also match gene boundaries
        first_exon_start = [exon.start for exon in db.children(gene, featuretype='exon', order_by='start')][0]
        last_exon_end    = [exon.end   for exon in db.children(gene, featuretype='exon', order_by='start')][-1]

        # check if gene boundaries correspond with CDS boundaries 
        first_cds_start = [cds.start for cds in db.children(gene, featuretype='CDS', order_by='start')][0]
        last_cds_end    = [cds.end   for cds in db.children(gene, featuretype='CDS', order_by='start')][-1]

        # if gene does have at least one UTR feature
        if 'three_prime_UTR' in gene_features or 'five_prime_UTR' in gene_features:

            # ## CDS boundaries should be within gene boundaries
            # if first_cds_start < gene.start:
            #     errors.append( ValueError(f"The start coordinate of {gene.id} on {gene.seqid} is not before the first CDS start coordinate") )
            # if last_cds_end > gene.end:
            #     errors.append( ValueError(f"The end coordinate of {gene.id} on {gene.seqid} is not after the last CDS end coordinate") )
        
            ## exon boundaries should be within gene boundaries
            if first_exon_start < gene.start:
                errors.append( ValueError(f"The start coordinate of {gene.id} on {gene.seqid} is not before the first exon start coordinate") )
            if last_exon_end > gene.end:
                errors.append( ValueError(f"The end coordinate of {gene.id} on {gene.seqid} is not after the last exon end coordinate") )

        # but if it does NOT have UTR features,
        else:

            # ## CDS boundaries should match gene boundaries
            # if first_cds_start != gene.start:
            #     errors.append( ValueError(f"The start coordinate of {gene.id} on {gene.seqid} does not correspond to the first CDS start coordinate") )
            # if last_cds_end != gene.end:
            #     errors.append( ValueError(f"The end coordinate of {gene.id} on {gene.seqid} does not correspond to the last CDS end coordinate") )
        
            ## exon boundaries should match gene boundaries
            if first_exon_start != gene.start:
                errors.append( ValueError(f"The start coordinate of {gene.id} on {gene.seqid} does not correspond to the first exon start coordinate") )
            if last_exon_end != gene.end:
                errors.append( ValueError(f"The end coordinate of {gene.id} on {gene.seqid} does not correspond to the last exon end coordinate") )


    if 'CDS' in gene_features:

        # check if gene boundaries correspond with CDS boundaries
        first_cds_start = [cds.start for cds in db.children(gene, featuretype='CDS', order_by='start')][0]
        last_cds_end    = [cds.end   for cds in db.children(gene, featuretype='CDS', order_by='start')][-1]

        if 'three_prime_UTR' in gene_features or 'five_prime_UTR' in gene_features:

            ## CDS boundaries should fall within gene boundaries, if we have UTR features
            if first_cds_start < gene.start:
                errors.append( ValueError(f"The start coordinate of {gene.id} on {gene.seqid} is not before the first CDS start coordinate") )
            if last_cds_end > gene.end:
                errors.append( ValueError(f"The end coordinate of {gene.id} on {gene.seqid} is not after the last CDS end coordinate") )

        else:
            # CDS boundaries should MATCH gene boundaries, if we don't have UTR features
            if first_cds_start != gene.start:
                errors.append( ValueError(f"The start coordinate of {gene.id} on {gene.seqid} does not correspond to the first CDS start coordinate") )
                continue
            if last_cds_end != gene.end:
                errors.append( ValueError(f"The end coordinate of {gene.id} on {gene.seqid} does not correspond to the last CDS end coordinate") )
                continue

        # store nt sequence in 'seq'
        seq = ''

        # if multiple CDS features,
        # how to merge them depends on strand

        # if + strand, the order is smallest-to-largest coordinate CDSs
        if gene.strand == '+':
            for cds in db.children(gene, featuretype='CDS', order_by='start'):
                seq += cds.sequence(fasta)

        # if - strand, the order is largest-to-smallest coordinate CDSs
        # NOTE: the nt sequence does not need to be revcomp'd,
        # this is probably taken care of by cds.sequence(fasta) ; the 'cds' object has strand information`
        elif gene.strand == '-':
            for cds in db.children(gene, featuretype='CDS', order_by='start', reverse=True):
                seq += cds.sequence(fasta)
                
            
        # the CDS nucleotide sequence should be a multiple of 3!
        if len(seq) % 3 != 0:
            errors.append( ValueError(f"The CDS nucleotide sequence length of {gene.id} on {gene.seqid} is not a multiple of 3! Remainder: {len(seq) % 3}") )

        # translate into protein CDS
        seq_record = Seq(seq)
        translation = str( seq_record.translate() )

        if translation[0] != 'M':
            errors.append( ValueError(f"The translation of {gene.id} on {gene.seqid} does not start with a regular START codon!") )

        if translation[-1] != '*':

            # if checking a blastocystis genome, we'd like to check for
            # the STOP codon motif T....TGTTTGTT if it does not have a regular STOP codon
            if args.is_blastocystis == True:

                # we need to extend gene coordinates to look for the motif
                if gene.strand == '+':
                    gene.end += 12 
                elif gene.strand == '-':
                    gene.start -= 12 

                # we only need to check the last bit of the gene for the motif
                motif_space_seq = gene.sequence(args.fasta_file)[-15:].upper()
                # print(f'{gene.id = }', file=sys.stderr)
                # print(f'{gene.start = }, {gene.end = }', file=sys.stderr)
                # print(f'{motif_space_seq = }', file=sys.stderr)

                # gradually allow for more substitutions
                for i in range(0,5):
                    # print(i)
                    # if re.match(r'(T|TA|TG)[ACTG]{4}((?:TGTTTGTT){s<=1})', motif_space_seq):
                    # pattern = re.compile(rf'(T|TA|TG)[ACTG]{{4}}((?:TGTTTGTT){{s<={i}}})')
                    pattern = re.compile(rf'(T|TA|TG)[ACTG]{{3,5}}((?:TGTTTGTT){{s<={i}}})')
                    # print(pattern.pattern)

                    # if 0 or 1 or 2 substitutions found a match, we believe its a genuine motif
                    if i <= 2 and pattern.search(motif_space_seq):
                        print(f'{gene.id} lacks a STOP codon but has a T....TGTTTGTT motif, distance of {i} substitutions', file=sys.stderr)
                        break
                    # if 3 substitutions found a match, maybe check it manually
                    elif i == 3 and pattern.search(motif_space_seq):
                        errors.append( ValueError(f'{gene.id} lacks a STOP codon but has a T....TGTTTGTT motif, distance of {i} substitutions') )
                        break
                    elif i == 4:
                        # if no matches were found even at 3 allowed substitutions, report no motif found
                        errors.append( ValueError(f"The translation of {gene.id} on {gene.seqid} does not end with a regular STOP codon and does not have a T....TGTTTGTT motif!") )

            # for any other organism, just report the missing STOP codon
            else:
                errors.append( ValueError(f"The translation of {gene.id} on {gene.seqid} does not end with a regular STOP codon") )

        # remove the trailing *
        cds_seq = translation.rstrip('*')

        # the CDS sequence should not have a * after stripping away the STOP codon's *
        if '*' in cds_seq:
            errors.append( ValueError(f"The CDS sequence of {gene.id} on {gene.seqid} contains at least one premature STOP codon!") )

        if args.print_proteins:
            print()
            print(f'>{gene.id}')
            print( cds_seq )

# print errors if we have them
if errors:
    for e in errors:
        print(e)
else:
    print("Gene models in provided GFF3 passed all checks!", file=sys.stderr)
