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

# standard stuff
import sys
import argparse

# stuff you need to install
import pysam
import regex as re
import gffutils

'''
Install required packages with
$ conda create -n <env_name> pysam regex gffutils
'''

# argument parser
parser = argparse.ArgumentParser(
            description=
        '''
        Find polyA sites that are supported by RNAseq data.

        PolyA tails in mRNA molecules are also sequenced in a
        typical RNAseq experiment and should be present in
        the reads.

        Obviously, polyA/T tails do not map to the genome, and hence
        these are found in the 'soft-clipped' part of a read alignment

        In reads from forward transcripts, polyA is at the end of the read mapping
        In reads from reverse transcripts, polyT is at the start of the read mapping

        This script looks for polyA/T tails of at least length 5
        in the soft-clipped parts of read alignments,
        and reports the last non-clipped part of the alignment as
        a polyA site as a GFF3 feature
        ''',
        formatter_class=argparse.RawTextHelpFormatter
        )

parser.add_argument(
        "-b", "--bam",
        dest='bam_file',
        type=str,
        required=True,
        help="Input RNAseq BAM file")

args = parser.parse_args()

# bam_file = sys.argv[1]

# r1_adapter: str = 'AGATCGGAAGAGCGTCGTGTAG'
# r2_adapter: str = 'AGATCGGAAGAGCACACGTCTG'

def main() -> None:
    # load BAM file
    # bam = pysam.AlignmentFile(bam_file, "rb")
    bam = pysam.AlignmentFile(args.bam_file, "rb")

    # container for found polyA sites
    found_sites = set()

    for contig in bam.references:
    # for contig in ['tig00000498']:

        print(f'Parsing {contig = } ...', file=sys.stderr)

        n = 0

        # for read in bam.fetch(contig, 1, 100000):
        for read in bam.fetch(contig):
            
            if not read.has_tag('XS'):
                continue

            if read.get_tag('XS') == '+':

                # read.cigartuples[-1][0] returns the last operation
                # (match, ins, del, softclip, etc) of the CIGAR string
                # 4 codes for Soft Clipping
                operation: int        = read.cigartuples[-1][0]
                soft_clip_length: int = read.cigartuples[-1][1]

                # skip if unmapped, not soft=clipped or short clip
                if read.is_unmapped or operation != 4 or soft_clip_length < 5:
                    continue

                # get the clipped sequence (clipped is at the end for '+' transcript reads)
                sequence: str = read.query_sequence[-soft_clip_length:].upper()

                # # skip if clipped sequence matches illumina adapter
                # pattern = re.compile(rf'(?:{sequence}){{e<=2}}')
                # if pattern.search(r2_adapter):
                #     continue

                # skip if sequence has anything but A's
                if re.search(r'[^A]', sequence):
                    continue

                found_sites.add( ('+', read.reference_end) )


            elif read.get_tag('XS') == '-':

                # read.cigartuples[0][0] returns the first operation
                # (match, ins, del, softclip, etc) of the CIGAR string
                # 4 codes for Soft Clipping
                operation: int        = read.cigartuples[0][0]
                soft_clip_length: int = read.cigartuples[0][1]

                # skip if unmapped, not soft=clipped or short clip
                if read.is_unmapped or operation != 4 or soft_clip_length < 5:
                    continue

                # get the clipped sequence (clipped is at the start for '-' transcript reads)
                sequence: str = read.query_sequence[:soft_clip_length].upper()

                # skip if sequence has anything but T's
                if re.search(r'[^T]', sequence):
                    continue

                found_sites.add( ('-', read.reference_start) )

        # generate GFF3 features for found polyA sites
        for site in found_sites:
            strand, coordinate = site

            n += 1
            polyA_site = gffutils.feature.Feature(
                    seqid = contig,
                    source = 'find_polyA_sites.py',
                    featuretype = 'polyA',
                    start = coordinate if strand == '+' else coordinate+1,
                    end = coordinate if strand == '+' else coordinate+1,
                    score = '.',
                    strand = strand,
                    frame = '.',
                    attributes = f'ID=polyA_site_{n:05d}_{strand}'
                    )
            print(polyA_site)


if __name__ == '__main__':
    main()

