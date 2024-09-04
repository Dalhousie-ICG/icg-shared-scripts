#!/usr/bin/env python

"""
Copyright 2023 Joran Martijn & Greg Seaton

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

import gffutils
import pysam
from statistics import mean
import re
from Bio import SeqIO
import random
import orfipy_core
import argparse
import os, subprocess
# from tqdm import tqdm


# argument parser
parser = argparse.ArgumentParser(
            description=
        '''
        Identifies regions between annotated genes and
        tries to predict genes there.

        Introns are predicted with bam2hints and
        filterIntronsFindStrand.pl from the Braker2 pipeline.

        Only introns that overlap "intergenic" regions and
        are well supported by the RNAseq data, as well fall
        within a certain user-specified size window, are considered.
        '''
        )
parser.add_argument(
        "-g", "--gff3",
        dest='gff3_file',
        type=str,
        required=True,
        help="Input GFF3 genome file")
parser.add_argument(
        "-b", "--bam",
        dest='bam_file',
        type=str,
        required=True,
        help="Input RNAseq BAM file")
parser.add_argument(
        "-f", "--fasta",
        dest='fasta_file',
        type=str,
        required=True,
        help="Input FASTA file")
parser.add_argument(
        "-w", "--window_size",
        dest='wdw_size',
        type=int,
        required=False,
        default=3,
        help="Window size for determining high-cov regions")
parser.add_argument(
        "-x", "--min_intron_length",
        dest='min_intron_len',
        type=int,
        required=False,
        default=10,
        help="Do not consider introns smaller than this")
parser.add_argument(
        "-y", "--max_intron_length",
        dest='max_intron_len',
        type=int,
        required=False,
        default=200,
        help="Do not consider introns larger than this")
parser.add_argument(
        "-l", "--min_orf_length",
        dest='min_orf_len',
        type=int,
        required=False,
        default=300,
        help="Do not consider orfs shorter than this")
args = parser.parse_args()



def get_igr_introns(db, intron_db, igrs):
    '''
    Find the introns that are in the intergenic regions
    '''
    igr_introns = []
    hold = []
    all_contigs = list( db.seqids() )
    # per contig, compare introns with intergenic regions
    # and select those that are within intergenic regions
    # for contig in tqdm( all_contigs, desc='Parsing contigs for intergenic introns' ):
    for contig in all_contigs:
        # select igrs and introns in a particular contig
        igrs_c    = (x for x in igrs if x.seqid == contig)
        introns_c = intron_db.region(seqid=contig, start=1)
        # compare all-vs-all
        for r in igrs_c:
            # check if introns in memory match the current region
            # if so, add them to list
            if len(hold) > 0:
                n=0
                for h in hold:
                    if h.end < r.start:
                        n+=1
                    elif h.start >= r.start and h.end <= r.end:
                        igr_introns.append(h)
                        n+=1
                hold = hold[n:]
            for i in introns_c:
                # select intron if within igr boundaries
                if i.start >= r.start and i.end <= r.end:
                    igr_introns.append(i)
                # skip all other introns comparisons once
                # we are past the considered igr
                if i.start > r.end:
                    # put intron that exceeded region
                    # into memory
                    hold.append(i)
                    break
                # note also that 'introns_c' is a generator object,
                # so each intron is only considered in a comparison once
                # this prevents many unnecessary comparisons
    return igr_introns


def get_true_introns(igr_introns, bamfile):
    '''
    Find the introns that are well supported by RNAseq data
    '''
    true_introns = []
    # loop over intergenic_region introns
    # for i in tqdm( igr_introns, desc='Getting well supported introns' ):
    for i in igr_introns:
        # start list of spliced / all aligned reads ratio's
        ratios = []
        # i.start and i.end are 1-indexed
        # pysam works with 0-indexing
        # so to get coverage for position x (1-indexed),
        # we need to ask for position x-1 (0-indexed)
        for col in bamfile.pileup(contig=i.seqid, start=i.start-1, stop=i.end, truncate=True, ignore_orphans=False):
            # skip col if no reads are aligned
            if col.get_num_aligned() == 0:
                continue
            # start spliced alignment counter
            spliced = 0
            # iterate over aligned reads at col
            for read in col.pileups:
                # add 1 to spliced alignment counter if read has 'N' here
                # in the CIGAR string. This indicates the read is spliced here
                if read.is_refskip:
                    spliced = spliced + 1
            # calculate ratio (spliced / all aligned reads)
            ratio = spliced / col.get_num_aligned()
            # append ratio to list of ratios
            ratios.append(ratio)

        # if ratio > 0.80,
        # intron is most likely true
        if len(ratios) > 0 and mean(ratios) > 0.50:
            i.attributes['ratio'] = f'{mean(ratios):.2f}'
            true_introns.append(i)
    return true_introns


def get_high_cov_regions(igr, bamfile, wdw_size):
    '''
    Find high coverage regions on the positive strand and
    the negative strand within the intergenic regions
    '''
    last_high_cov_pos = None
    for s in ['+', '-']:

        switch = False
        for col in bamfile.pileup(
                contig=igr.seqid,
                start=igr.start+20,
                stop=igr.end-20,
                truncate=True):

            # pysam 0-indexed, so add +1 to make it 1-indexed
            pos = col.reference_pos + 1
            # calculate site coverage, but stop iterating
            # once you have a coverage of 10
            # this saves a lot of computations
            cov = 0
            for x in col.pileups:
                if s in x.alignment.get_tag('XS'):
                    cov += 1
                if cov > 10:
                    break
            # cov = len( [x for x in col.pileups if s in x.alignment.get_tag('XS')] )
            # print(s, pos, ':', cov)
            # cov = col.nsegments
            # cov = get_strand_cov(col, strand)

            # if you go from low-cov to high-cov region
            if switch == False and cov >= 2:
                region_start = last_high_cov_pos = pos
                switch = True

            # if you are in a high-cov region and ...
            if switch == True:
                # ... encounter another high-cov position
                if cov >= 2:
                    # it may be after jumping over a long non-covered region
                    # in that case, yield the region prior to the jump
                    if pos - last_high_cov_pos > wdw_size:
                        region_end = last_high_cov_pos
                        yield (igr.seqid, region_start-25, region_end+25, s, region_end - region_start)
                        region_start = last_high_cov_pos = pos
                    # it may be the end of the intergenic region
                    # in that case, yield the current region
                    elif pos == igr.end-20:
                        region_end = pos
                        yield (igr.seqid, region_start-25, region_end+25, s, region_end - region_start)
                    # if it is a high-cov position close to the last one,
                    # reset the last_high_cov_pos
                    else:
                        last_high_cov_pos = pos

                # ... encounter a low-cov position
                else:
                    # and current pos is far away from the last highcov pos
                    # OR the current position is at the end of the igr,
                    # yield the last region before encountering low cov sites
                    if pos - last_high_cov_pos > wdw_size-1 or pos == igr.end-20:
                        region_end = last_high_cov_pos
                        yield (igr.seqid, region_start-25, region_end+25, s, region_end - region_start)
                        switch = False
                    # if current pos is close,
                    # keep switch on and do not set the region_end yet
                    else:
                        continue


def get_introns_of(region, introns):
    '''
    Find introns that are within a high coverage region
    '''
    # intron objects are 1-indexed
    # region is 1-indexed
    r_seqid, r_start, r_end, r_strand, r_length = region
    for i in introns:
        if i.seqid != r_seqid or i.strand != r_strand:
            continue
        if r_start < i.start and r_end > i.end:
            # transform coordinates so they are relative
            # to the high_cov_region start
            # 1-indexed - 1-indexed = 0-indexed
            # thus +1 to return to 1-indexing
            i.start = i.start - r_start + 1
            i.end   = i.end   - r_start + 1
            yield i
        elif i.start > r_end:
            break


def select_non_overlapping_introns(introns):
    '''
    Find the set of introns that do not overlap
    '''
    overlapping_introns = []
    prev_intron_start, prev_intron_end = 0, 100000000
    # iterate over introns one by one
    # assume they are sorted by start, then by end coords
    for i in introns:
        # if current intron overlaps with previous intron,
        # add it to overlapping_introns group
        if prev_intron_start <= i.start <= prev_intron_end:
            overlapping_introns.append(i)
        # if current intron does not overlap previous intron,
        # yield the best scored intron from the overlap group
        # and reset the intron group
        elif prev_intron_end < i.start:
            best_intron = max( overlapping_introns, key=lambda x : int(x.score) )
            yield best_intron
            overlapping_introns = [i]
        # update start and end coords of previous intron
        prev_intron_start, prev_intron_end = i.start, i.end
    # for the last intron group, select the best intron
    # after all introns have been iterated over
    if len(introns) > 0:
        best_intron = max( overlapping_introns, key=lambda x : x.score )
        yield best_intron


def splice_seq(seq, introns):
    '''
    Take a sequence string and an iterator or list of introns
    and splice out the introns from the sequence string
    '''

    # -1 from intron start (='exon' end)
    # +1 from intron end   (='exon' start)
    coords = [ [i.start-1 , i.end+1] for i in introns ]

    # flatten list and convert to 0-indexing
    f = [ x-1 for c in coords for x in c ]
    # add start and end (0-indexed)
    f.insert(0, 0)
    f.append(len(seq)-1)

    # repack list
    borders = [ [f[i],f[i+1]] for i in range(0, len(f), 2) ]
    # construct spliced sequence
    spliced_seq = ''.join([ seq[ b[0]:b[1]+1 ] for b in borders ])
    return spliced_seq


def get_non_overlapping_orfs(seq, strand):
    '''
    Given a DNA sequence string, predict the ORFs with orfipy
    and select the biggest ORFs that do not overlap

    Modified for Blastocystis to find STOP codons for genes
    with the T....TGTTTGTT motif
    '''
    seq = seq.upper()
    # convert strand symbol
    if strand == '+':
    #     # transform stop motifs into stop codons
    #     # (a Blastocyctis specific thing)
    #     seq = re.sub(r'([ACTG]*)T[ACTG]{2}([ACTG]{2}TGTTTGTT)', r'\1TAA\2', seq)
        strand = 'f'
    elif strand == '-':
    #     # transform stop motifs into stop codons (reverse direction)
    #     # (a Blastocyctis specific thing)
    #     seq = re.sub(r'(AACAAACA[ACTG]{2})[ACTG]{2}A', r'\1TTA', seq, count=1)
        strand = 'r'
    # find orfs in seq
    ## orfs() returns a list of tuples
    all_orfs = orfipy_core.orfs(
        seq,
        minlen=args.min_orf_len,
        maxlen=100000,
        starts=['ATG'],
        strand=strand,
        include_stop=True
    )
    # first sort orfs by length
    all_orfs.sort(key=lambda x : x[1]-x[0], reverse=True)
    # then select non-overlapping orfs,
    # preferring the largest ones
    occupied_regions = []
    for i, orf in enumerate(all_orfs):
        # get coordinates, frame and length of orf
        start = orf[0]
        stop  = orf[1]
        # always select the first orf,
        # i.e. the largest orf
        if i == 0:
            occupied_regions.append([start,stop])
            yield orf
            continue
        # check if current orf overlaps with
        # any of the already selected orfs
        for region in occupied_regions:
            # move to next orf if any region overlaps
            if max(region[0], start) < min(region[1], stop):
                break
        # this else block is only executed if the above
        # for loop is exited normally (i.e. doesn't break)
        else:
            # if it survived the overlap checks, select orf
            occupied_regions.append([start,stop])
            yield orf


def re_insert_introns_to(orfs, introns, region_start):
    '''
    Adjust start and end coordinates of predicted ORFs by
    checking which introns overlap with the ORFs
    '''
    # orfs is a list of tuples (0-indexed),
    # introns is a list of intron objects (1-indexed)
    # convert to lists of lists so I can edit values
    upd_orfs = [ list(o) for o in orfs ]
    for o in upd_orfs:
        for i in introns:
            # junction is the coord in the spliced seq
            # where the intron is supposed to be
            # it is updated with each intron we visit / re-insert
            junction = i.start - 1 - 1 # -1 to make 0-indexed, -1 to get to junction
            intron_len = len(i)
            # if junction is situated before predicted orf
            if junction < o[0]:
                # update orf start and end coords
                o[0] += intron_len
                o[1] += intron_len
            # if junction is within the predicted orf
            elif o[0] <= junction < o[1]-1:
                # update only the orf end coord
                o[1] += intron_len
            # if junction is after the predicted orf
            elif o[1]-1 < junction:
                # move on to the next orf
                break
        # transform back to genome space
        o[0] = o[0] + region_start
        o[1] = o[1] + region_start - 1
    return upd_orfs


def make_features(region, orfs, introns, new_orf_count):
    '''
    Given a set of new ORFs in a particular region,
    and a set of overlapping introns, create
    new gene and CDS Feature objects
    '''
    possible_phases = [0,2,1]
    r_seqid, r_start, r_end, r_strand, _ = region

    # convert intron coordinates to genome space
    for i in introns:
       i.start = i.start + r_start - 1
       i.end   = i.end   + r_start - 1

    # move through orfs-vs-introns to make features
    for o in orfs:
        gene_start  = o[0]
        gene_end    = o[1]
        gene_strand = o[2]
        # initialize exon_start,
        # it will be updated as we move through introns
        exon_start = gene_start
        exon_end   = gene_end
        # create new gene id
        # new_gene_id = create_new_id(new_feature_ids)
        # new_feature_ids.append(new_gene_id)
        new_gene_name = f'orf_added{new_orf_count:04d}'
        contig_number = int( r_seqid.split('tig')[1] )
        new_gene_contig = f'ctg{contig_number:03d}'
        # create Feature obj
        new_gene = gffutils.feature.Feature(
                seqid = r_seqid,
                source = 'script',
                featuretype = 'gene',
                start = gene_start,
                end = gene_end,
                score = '.',
                strand = gene_strand,
                frame = '.',
                attributes = f'ID={new_gene_contig}.{new_gene_name};Name={new_gene_name}'
                )
        yield new_gene

        # the mRNA has the same coords as the gene
        new_mRNA_name = new_gene_name.replace('orf','mRNA')
        new_mRNA = gffutils.feature.Feature(
                seqid = r_seqid,
                source = 'script',
                featuretype = 'mRNA',
                start = gene_start,
                end = gene_end,
                score = '.',
                strand = gene_strand,
                frame = '.',
                attributes = f'ID={new_gene_contig}.{new_mRNA_name};Parent={new_gene_contig}.{new_gene_name};Name={new_mRNA_name}'
                )
        yield new_mRNA

        # if the region has no introns
        if len(introns) == 0:

            # the exon also has the same coords as the gene
            new_exon = gffutils.feature.Feature(
                    seqid = r_seqid,
                    source = 'script',
                    featuretype = 'exon',
                    start = gene_start,
                    end = gene_end,
                    score = '.',
                    strand = gene_strand,
                    frame = '.',
                    attributes = f'ID={new_gene_contig}.{new_mRNA_name}.exon01;Parent={new_gene_contig}.{new_mRNA_name}'
                    )
            yield new_exon

            # the CDS has the same coords as the gene
            new_cds = gffutils.feature.Feature(
                    seqid = r_seqid,
                    source = 'script',
                    featuretype = 'CDS',
                    start = gene_start,
                    end = gene_end,
                    score = '.',
                    strand = gene_strand,
                    frame = str(0),
                    attributes = f'ID={new_gene_contig}.{new_mRNA_name}.CDS01;Parent={new_gene_contig}.{new_mRNA_name}'
                    )
            yield new_cds


        # else if new gene does have introns,
        # CDS feature coords need to be adjusted
        else:
            # phase of first CDS is always 0
            # but adjustment of phase depends on
            # the strand of the gene
            phase = 0
            m = 1 # exon counter
            if gene_strand == '+':
                for i in introns:
                    if i.end < gene_start:
                        continue
                    elif gene_start < i.end < gene_end:
                        cds_start = exon_start
                        cds_end   = i.start - 1
                        new_exon = gffutils.feature.Feature(
                                seqid = r_seqid,
                                source = 'script',
                                featuretype = 'exon',
                                start = cds_start,
                                end = cds_end,
                                score = '.',
                                strand = '+',
                                frame = '.',
                                attributes = f'ID={new_gene_contig}.{new_mRNA_name}.exon{m:02d};Parent={new_gene_contig}.{new_mRNA_name}'
                                )
                        yield new_exon
                        new_cds = gffutils.feature.Feature(
                                seqid = r_seqid,
                                source = 'script',
                                featuretype = 'CDS',
                                start = cds_start,
                                end = cds_end,
                                score = '.',
                                strand = '+',
                                frame = str(phase),
                                attributes = f'ID={new_gene_contig}.{new_mRNA_name}.CDS{m:02d};Parent={new_gene_contig}.{new_mRNA_name}'
                                )
                        yield new_cds
                        new_intron = gffutils.feature.Feature(
                                seqid = r_seqid,
                                source = 'script',
                                featuretype = 'intron',
                                start = i.start,
                                end = i.end,
                                score = '.',
                                strand = r_strand,
                                frame = '.',
                                attributes = f'ID={new_gene_contig}.{new_mRNA_name}.intron{m:02d};Parent={new_gene_contig}.{new_mRNA_name}'
                                )
                        yield new_intron

                        # update exon start and phase for the next exon/cds
                        exon_start = i.end + 1
                        codon_remainder = ( len(new_cds) - phase ) % 3
                        phase = possible_phases[codon_remainder]
                        m += 1
                    elif gene_end < i.start:
                        break

                # update the last exon/cds coords
                cds_start = exon_start
                cds_end   = gene_end
                new_exon = gffutils.feature.Feature(
                        seqid = r_seqid,
                        source = 'script',
                        featuretype = 'exon',
                        start = cds_start,
                        end = cds_end,
                        score = '.',
                        strand = '+',
                        frame = '.',
                        attributes = f'ID={new_gene_contig}.{new_mRNA_name}.exon{m:02d};Parent={new_gene_contig}.{new_mRNA_name}'
                        )
                yield new_exon
                new_cds = gffutils.feature.Feature(
                        seqid = r_seqid,
                        source = 'script',
                        featuretype = 'CDS',
                        start = cds_start,
                        end = cds_end,
                        score = '.',
                        strand = '+',
                        frame = str(phase),
                        attributes = f'ID={new_gene_contig}.{new_mRNA_name}.CDS{m:02d};Parent={new_gene_contig}.{new_mRNA_name}'
                        )
                yield new_cds

            elif gene_strand == '-':
                # introns.reverse()
                for i in introns[::-1]:
                    if gene_end < i.start:
                        continue
                    elif gene_start < i.end < gene_end:
                        cds_start = i.end + 1
                        cds_end   = exon_end
                        new_exon = gffutils.feature.Feature(
                                seqid = r_seqid,
                                source = 'script',
                                featuretype = 'exon',
                                start = cds_start,
                                end = cds_end,
                                score = '.',
                                strand = '-',
                                frame = '.',
                                attributes = f'ID={new_gene_contig}.{new_mRNA_name}.exon{m:02d};Parent={new_gene_contig}.{new_mRNA_name}'
                                )
                        yield new_exon
                        new_cds = gffutils.feature.Feature(
                                seqid = r_seqid,
                                source = 'script',
                                featuretype = 'CDS',
                                start = cds_start,
                                end = cds_end,
                                score = '.',
                                strand = '-',
                                frame = str(phase),
                                attributes = f'ID={new_gene_contig}.{new_mRNA_name}.CDS{m:02d};Parent={new_gene_contig}.{new_mRNA_name}'
                                )
                        yield new_cds

                        new_intron = gffutils.feature.Feature(
                                seqid = r_seqid,
                                source = 'script',
                                featuretype = 'intron',
                                start = i.start,
                                end = i.end,
                                score = '.',
                                strand = r_strand,
                                frame = '.',
                                attributes = f'ID={new_gene_contig}.{new_mRNA_name}.intron{m:02d};Parent={new_gene_contig}.{new_mRNA_name}'
                                )
                        yield new_intron

                        # update exon_end and phase for next exon/cds
                        exon_end = i.start - 1
                        codon_remainder = ( len(new_cds) - phase ) % 3
                        phase = possible_phases[codon_remainder]
                        m += 1
                    elif i.end < gene_start:
                        # introns.reverse()
                        break

                # update the last exon/cds coords
                cds_start = gene_start
                cds_end   = exon_end
                new_exon = gffutils.feature.Feature(
                        seqid = r_seqid,
                        source = 'script',
                        featuretype = 'exon',
                        start = cds_start,
                        end = cds_end,
                        score = '.',
                        strand = '-',
                        frame = '.',
                        attributes = f'ID={new_gene_contig}.{new_mRNA_name}.exon{m:02d};Parent={new_gene_contig}.{new_mRNA_name}'
                        )
                yield new_exon
                new_cds = gffutils.feature.Feature(
                        seqid = r_seqid,
                        source = 'script',
                        featuretype = 'CDS',
                        start = cds_start,
                        end = cds_end,
                        score = '.',
                        strand = '-',
                        frame = str(phase),
                        attributes = f'ID={new_gene_contig}.{new_mRNA_name}.CDS{m:02d};Parent={new_gene_contig}.{new_mRNA_name}'
                        )
                yield new_cds

        # update orf count for each new orf in the igr
        new_orf_count += 1


# run bam2hints and filterIntronsFindStrand.pl if their output doesn't exist already
if not os.path.exists('scored_intron_hints.gff'):
    # print("Generating hints file")
    # run bam2hints
    subprocess.run(["bam2hints", "--intronsonly", f"--minintronlen={args.min_intron_len}", f"--in={args.bam_file}", "--out=intron_hints.gff"])
    # run filterIntronsFindStrand.pl and select introns with score 5 or higher
    cmd = "filterIntronsFindStrand.pl " + args.fasta_file + " intron_hints.gff --score --allowed=gtag,gcag,atac,atag,gaag | awk \'$6 >= 5 {print $0}\' > scored_intron_hints.gff"
    # cmd = f'filterIntronsFindStrand.pl {args.fasta_file} intron_hints.gff --score --allowed=gtag,gcag,atac,atag,gaag | awk "$6 >= 5 {print $0}" > scored_intron_hints.gff'
    subprocess.run([cmd], shell=True)
    # print("Hints file complete")

# load gff file
## intergenic regions do not necessarily have to be explicitly defined
## this script will infer them if they are not defined already
db = gffutils.create_db(args.gff3_file, ':memory:', merge_strategy='create_unique')

# load intron gff file
intron_db = gffutils.create_db('scored_intron_hints.gff', ':memory:')

# load bam file
bamfile = pysam.AlignmentFile(args.bam_file, mode='rb')

# get all gene features
genes = db.features_of_type('gene', order_by=('seqid','start') )
# infer intergenic_region features
igrs  = list( db.interfeatures(genes, new_featuretype='intergenic_region') )

# get introns in intergenic regions
igr_introns = get_igr_introns(db, intron_db, igrs)
# select introns of a certain size
sized_introns = list( filter(lambda x : args.min_intron_len < len(x) < args.max_intron_len, igr_introns) )
# have at least an x amount of support
true_introns = get_true_introns(sized_introns, bamfile)


# load fasta
fasta = SeqIO.index(args.fasta_file, 'fasta')

# main code
# new_feature_id_list = []
new_gene_count = 1
for igr in igrs:
    if len(igr) < 500:
        continue

    contig = fasta[igr.seqid].seq

    # print()
    # print(f'Intergenic region {igr.seqid} : {igr.start} - {igr.end}')

    igr_high_cov_regions = get_high_cov_regions(igr, bamfile, args.wdw_size)
    for r in igr_high_cov_regions:
        r_seqid, r_start, r_stop, r_strand, _ = r
        if r_stop - r_start < 300:
            continue
        # print()
        # print(f'high_cov_region: {r}')
        r_introns = list( get_introns_of(r, true_introns) )
        r_n_introns = list( select_non_overlapping_introns(r_introns) )
        r_seq = str( contig[ r_start-1:r_stop ] )
        r_spliced_seq = splice_seq(r_seq, r_n_introns)
        # print(r_spliced_seq)
        r_orfs = list( get_non_overlapping_orfs(r_spliced_seq, r_strand) )
        # print(r_orfs)
        if len(r_orfs) == 0:
            continue
        r_upd_orfs = re_insert_introns_to(r_orfs, r_n_introns, r_start)
        # print(r_upd_orfs)
        # r_new_features = make_features(r, r_upd_orfs, r_n_introns, new_feature_id_list)
        r_new_features = list( make_features(r, r_upd_orfs, r_n_introns, new_gene_count) )
        new_gene_count += len( [f for f in r_new_features if f.featuretype=='gene'] ) 
        db.update(r_new_features, merge_strategy='create_unique')

# print(db.directives)
# print directives/pragmas
for p in db.directives:
    print(f'##{p}')

# print some logging info
print(
    f'#This GFF3 file was created with {__file__}\n'
    f'#--bam {args.bam_file}'
    f' --gff {args.gff3_file}'
    f' --window_size {args.wdw_size}'
    f' --min_intron_length {args.min_intron_len}'
    f' --max_intron_length {args.max_intron_len}'
    f' --min_orf_length {args.min_orf_len}'
)

# print updated gff record in logical order
for g in db.features_of_type('gene', order_by=('seqid', 'start')):
    print()
    print(g)
    for f in db.children(g, order_by='start'):
        print(f)

