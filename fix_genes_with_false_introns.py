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
from statistics import mean
import re
import sys

# stuff you need to install
import gffutils
import pysam
from Bio.Seq import Seq
from Bio import SeqIO
import orfipy_core


# make command line interface
parser = argparse.ArgumentParser(
        description=
        """
        IMPORTANT: This script expects that there are explicit intron
        features in the GFF3 file!! If you're GFF3 doesn't have them,
        consider adding them with add_intron_features.py!!

        Checks whether predicted introns are supported by RNAseq data,
        and removes them if they are not.

        If a gene is found to have an intron not supported by the
        RNAseq data, it and all its daughter features (including
        all supported introns) are stripped away.

        The script then attempts to define a region in which to find
        new ORFs. 
            By default, the region starts at the end of the 
        previous gene +1, and ends at the start of the next gene -1.
            If the deleted gene overlapped with previous or next gene,
        the region will start & end at the end/start of the first 
        non-overlapping upstream/downstream gene, respectively.
            If it can't find a non-overlapping upstream/downstream gene,
        (for example because it is close to the contig edge), the
        region will start at the start of the most upstream gene,
        or end at the end of the most downstream gene.
            If the deleted gene is the first gene of the contig, it will
        start at the start of the gene and end at the start of the next gene -1
        If the deleted gene is the last gene of the contig, it will
        start at the end of the previous gene +1 and end at the end
        of the deleted gene.
            If the deleted gene directly followed another deleted gene,
        the script will keep track of any new ORFs that were inferred
        while it was visiting that other deleted gene, and adjust the
        start of the region accordingly.

        The script then attempts to re-evaluate the intron locations
        from the RNAseq BAM file in the defined region. If any two
        introns overlap, it will select the intron with the most support.
        The introns are spliced out,
        (one spliced sequence per strand), and predict a set of
        non-overlapping ORFs on each strand. Both sets of ORFs are then
        pooled, and the biggest ORFs are selected, allowing for partial
        overlap between ORFs of different strands.
        
        Any true intron that overlaps with the newly selected ORFs are
        then re-inserted, and new gene, mRNA, CDS, exon and intron
        features are created and written to a new GFF3 file.
        """,
        formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument(
    "-g", "--gff3",
    dest="gff3_file",
    metavar="GFF3_FILE",
    required=True,
    help="Input GFF3 file. Must contain explicitly written intron features")

parser.add_argument(
    "-b", "--bam_file",
    dest="bam_file",
    metavar="BAM_FILE",
    required=True,
    help="Input BAM file with RNAseq mapping")

parser.add_argument(
    "-f", "--fasta_file",
    dest="fasta_file",
    metavar="FASTA_FILE",
    required=True,
    help="Input FASTA file with genome sequence")

parser.add_argument(
    "-t", "--threshold",
    dest="threshold",
    metavar="THRESHOLD",
    type=float,
    default=0.50,
    required=False,
    help="Intron support threshold. Value between 0 and 1. Default 0.50")

parser.add_argument(
    "-s", "--splice_sites",
    dest="splice_sites",
    nargs='+',
    type=str,
    default=["gtag", "gcag", "atac", "atag", "gaag", "acac", "ggag"],
    required=False,
    help="Allow only introns that have these types of splice sites")


args = parser.parse_args()


def main(args) -> None:

    # load gff file with introns explicitly defined
    db = gffutils.create_db(args.gff3_file, ':memory:', merge_strategy='create_unique')

    # load bam file
    bamfile = pysam.AlignmentFile(args.bam_file, mode='rb')

    # load fasta
    fasta = SeqIO.index(args.fasta_file, 'fasta')

    # main code
    false_genes = []
    new_features = []

    # treat each contig seperately
    for contig in db.seqids():

        # reset furthest gene end
        furthest_gene_end = 0

        # put all gene objects in memory,
        # we need this to look up previous and next gene starts/ends
        contig_genes = list( db.region(seqid=contig, start=1, featuretype='gene') )

        for i, gene in enumerate( contig_genes ):#, start=1 ):
            # if '000090' not in gene.id: continue
            # print()
            # print(i, gene)

            # get introns belonging to the gene
            introns = db.children(gene, featuretype='intron', order_by='start')

            # introns are 1-indexed
            # get_intron_rnaseq_support() expects 0-indexed
            supports = [ get_intron_rnaseq_support(j.seqid, j.start-1, j.end-1, bamfile) for j in introns ]

            # if all introns are supported, do nothing and move to next gene
            if all( support >= args.threshold for support in supports ): continue

            # mark gene for deletion
            false_genes.append(gene)

            # delete all introns related to this gene
            db.delete(introns)

            # region = previous gene end + 1 and next gene start - 1
            # get region coords and region sequence

            ## if a new orf was predicted beyond the original end of the previous gene,
            ## 'furthest_gene_end' will have a larger value than previous_gene_end + 1
            ## region_start will thus be furthest_gene_end + 1 instead
            # print(f'{furthest_gene_end = }')

            region_start, region_end = get_region(contig_genes, i, furthest_gene_end)
            # print('region: ', region_start, region_end)

            # look for new genes on both strands independently
            new_orfs = []

            # new_introns = []
            new_introns = {}

            # if an intron on one strand overlaps with the region_start,
            # the region_start will be different for each strand
            # so we need check for that and agree on a region_start before
            # predicting ORFs etc
            for strand in ['+', '-']:
                # get introns directly from bam for this region
                ## region_start and region_end are 1-indexed
                ## bamfile.fetch() expects 0-indexed
                selected_introns = get_introns_from_bam(bamfile, contig, strand, region_start-1, region_end-1, fasta)
                # print()
                # print('selected_introns')
                # for intron in selected_introns:
                #     print(intron)

                # ensure introns do not overlap
                sorted_introns = sorted( selected_introns, key=lambda x: x.start )
                disjoint_introns = list( select_non_overlapping_introns(sorted_introns) )
                # print()
                # print('disjoint_introns')
                # for intron in disjoint_introns:
                #     print(intron)

                # adjust region_start and/or region_end if
                # they are inside an intron, such that the 
                # partially overlapping intron is excluded
                for i in disjoint_introns:
                    if i.start < region_start < i.end:
                        region_start = i.end+1
                        # disjoint_introns.remove(i)
                    elif i.start < region_end < i.end:
                        region_end   = i.start-1
                        # disjoint_introns.remove(i)
                disjoint_introns = [ j for j in disjoint_introns if region_start < j.start < j.end < region_end ]
                new_introns[strand] = disjoint_introns
                # print(f'{region_start = }')

            # transform intron coords to region_seq space 
            # for i in disjoint_introns:
            for strand in ['+','-']:
                # print(f'{strand = }')
                for i in new_introns[strand]:
                    # print(f'{i = }')
                    i.start = i.start-region_start+1 
                    i.end   = i.end  -region_start+1 
                    # print(f'{region_start = }')

            # splice out introns and predict ORFs
            for strand in ['+','-']:

                # retrieve contig DNA sequence
                contig_seq = fasta[contig].seq
                ## python slicing is 0-indexed, so substract 1
                ## from 1-indexed region_start value
                region_seq = str( contig_seq[region_start-1:region_end] )

                # splice introns from region_seq
                # spliced_region_seq  = splice_seq(region_seq, disjoint_introns)
                spliced_region_seq  = splice_seq(region_seq, new_introns[strand])

                # predict orfs, only on strand of original gene
                strand_new_orfs = get_non_overlapping_orfs(spliced_region_seq, strand)

                # adjust orf start and stop by re-inserting overlapping introns
                # strand_intron_corrected_orfs = re_insert_introns(strand_new_orfs, disjoint_introns)
                strand_intron_corrected_orfs = re_insert_introns(strand_new_orfs, new_introns[strand])
                # for corrected_orf in strand_intron_corrected_orfs:
                #     print(f'{corrected_orf = }')

                # add identified disjoint introns to new_introns list
                # new_introns.extend(disjoint_introns)
                # add identified orfs to new_orfs list
                new_orfs.extend(strand_intron_corrected_orfs)


            # skip to next gene if no orfs were found
            if len(new_orfs) == 0: continue

            # select biggest ORFs, and allow for partial overlap between genes of opposing strands
            selected_orfs = select_orfs(orfs=new_orfs, overlap='partial')
            # print()
            # print(f'{selected_orfs = }')
            
            # transform orfs back to genome space
            region_orfs = [ [o[0]+region_start,
                             o[1]+region_start-1,
                             o[2],o[3]] 
                           for o in selected_orfs ]
            # sort by starting coordinate
            region_orfs.sort(key=lambda x: x[0])
            # print()
            # for orf in region_orfs:
            #     print('region orf: ', orf)


            # transform intron coords back to genome space 
            for _ in new_introns.values():
                for i in _:
                    i.start = i.start+region_start-1 
                    i.end   = i.end  +region_start-1 
                    # print(f'{i = }')

            # loop over orfs and create new gene and mRNA for each new orf
            n = 1 # gene counter
            for orf in region_orfs:

                # create new gene
                new_gene = create_new_feature('gene', orf, templ=gene, gene_number=n)

                # create new mRNA
                new_mrna = create_new_feature('mRNA', orf, templ=gene, gene_number=n)

                # add to new_features
                new_features.extend( [new_gene,new_mrna] )

                # update furthest gene end
                furthest_gene_end = max( new_gene.end, furthest_gene_end )

                # loop over introns and create new exons, CDSs and introns along the way
                phase = 0
                m = 1 # exon counter

                # if orf is on + strand
                if orf[2] == '+':
                    moving_start = orf[0]
                    # for i in new_introns:
                    for i in new_introns['+']:
                        # skip if intron doesn't match the strand
                        # if i.strand != '+':
                        #     continue
                        # print(f'{i = }')
                        # if intron before gene
                        if i.end < orf[0]:
                            # delete intron from database and move to next intron
                            db.delete(i)
                            continue
                        # if intron within gene
                        elif orf[0] < i.end < orf[1]:
                            # create new exon
                            exon_start  = moving_start
                            exon_end    = i.start - 1
                            new_exon = create_new_feature('exon', orf, templ=gene, gene_number=n, start=exon_start, end=exon_end, exon_number=m)
                            # print(f'{new_exon = }')
                            # create new CDS
                            new_cds = create_new_feature('CDS', orf, templ=gene, gene_number=n, start=exon_start, end=exon_end, phase=phase)

                            # add to new_features
                            new_features.extend( [new_exon,new_cds] )

                            # update phase, moving start and counter for next exon
                            phase = determine_new_phase(new_exon, phase)
                            moving_start = i.end + 1
                            m += 1
                        # if intron after gene
                        elif orf[1] < i.start:
                            # move to next orfs
                            break

                    # the last exon and CDS
                    # or first exon and CDS if there were no introns
                    exon_start  = moving_start
                    exon_end    = orf[1]
                    new_exon = create_new_feature('exon', orf, templ=gene, gene_number=n, start=exon_start, end=exon_end, exon_number=m)
                    # print(f'{new_exon = }')
                    new_cds = create_new_feature('CDS', orf, templ=gene, gene_number=n, start=exon_start, end=exon_end, phase=phase)

                    # add to new_features
                    new_features.extend( [new_exon,new_cds] )


                elif orf[2] == '-':
                    moving_end = orf[1]
                    # for i in new_introns[::-1]:
                    for i in new_introns['-'][::-1]:
                        # skip if intron doesn't match the strand
                        # if i.strand != '-':
                        #     continue
                        # print(f'{i = }')
                        # if intron after gene
                        if orf[1] < i.start:
                            continue
                        # if intron within gene
                        elif orf[0] < i.end < orf[1]:
                            # create new exon
                            exon_start  = i.end + 1
                            exon_end    = moving_end
                            new_exon = create_new_feature('exon', orf, templ=gene, gene_number=n, start=exon_start, end=exon_end, exon_number=m)
                            # print(f'{new_exon = }')
                            # create new CDS
                            new_cds = create_new_feature('CDS', orf, templ=gene, gene_number=n, start=exon_start, end=exon_end, phase=phase)

                            # add to new_features
                            new_features.extend( [new_exon,new_cds] )

                            # update phase, moving end and counter for next exon
                            phase = determine_new_phase(new_exon, phase)
                            moving_end = i.start - 1
                            m += 1
                        # if intron before gene
                        elif i.end < orf[0]:
                            break
                    # the last exon and CDS
                    # or first exon and CDS if there were no introns
                    exon_start  = orf[0] 
                    exon_end    = moving_end
                    new_exon = create_new_feature('exon', orf, templ=gene, gene_number=n, start=exon_start, end=exon_end, exon_number=m)
                    # print(f'{new_exon = }')
                    new_cds = create_new_feature('CDS', orf, templ=gene, gene_number=n, start=exon_start, end=exon_end, phase=phase)

                    # add to new_features
                    new_features.extend( [new_exon,new_cds] )

                # update gene counter
                n += 1


    # delete the old gene with false introns
    # and delete all its children
    for gene in false_genes:
        for f in db.children(gene, featuretype=('exon','intron','CDS','mRNA')):
            db.delete(f)
        db.delete(gene)

    # update db with all new_features
    db.update(new_features, merge_strategy='create_unique')

    # print directives/pragmas
    for p in db.directives:
        print(f'##{p}')

    # print some logging info
    print(
        f'#This GFF3 file was created with {__file__}\n'
        f'#--bam_file {args.bam_file}'
        f' --gff3 {args.gff3_file}'
        f' --fasta {args.fasta_file}'
        f' --threshold {args.threshold}'
    )

    # print updated gff record in logical order
    for gene in db.features_of_type('gene', order_by=('seqid', 'start')):
        print()
        print(gene)
        for f in db.children(gene, order_by='start'):
            print(f)

    pass


def get_intron_rnaseq_support(contig, start, stop, bamfile):
    '''
    Reports fraction of mapped read positions that
    indicate a spliced read (i.e. an intron)
    '''
    # print(f'Checking intron {contig}:{start}-{stop} for RNAseq support', file=sys.stderr)
    ratios = []
    columns = bamfile.pileup(contig=contig, start=start, stop=stop, truncate=True, ignore_orphans=False)

    # may be there are no reads whatsoever mapping
    # to this intron. In that case support is 0
    if not any(columns):
        return 0.0

    for col in columns: 
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
    # very rarely the columns is neither empty nor iterated over
    # (I don't really understand it), but then ratios does not
    # get populated. So in that case, return support 0
    if len(ratios) == 0:
        return 0.0
    return mean(ratios)


def get_region(
        contig_genes: list[gffutils.feature.Feature], 
        i: int, 
        furthest_gene_end: int
    ) -> tuple[int, int]:
    '''
    Get the start and end coordinates of a sequence region
    The region stretches from 
    the first base after the previous gene until
    the first base before the next gene
    '''

    # calculate total number of genes on contig
    gene_number = len(contig_genes)

    # if first gene on contig,
    if i == 0:
        region_start = contig_genes[i].start
        region_end   = contig_genes[i+1].start-1

    # if last gene on contig,
    # region starts at first base after previous gene
    # region ends   at last  base of    current  gene
    elif i == gene_number-1:
        region_start = max( contig_genes[i-1].end+1, furthest_gene_end+1 )
        region_end   = contig_genes[i].end

    # if the current gene overlaps with the previous gene,
    elif contig_genes[i-1].start < contig_genes[i].start < contig_genes[i-1].end:

        # iterate over upstream genes
        for j in range(1, i+1):
            # if the most upstream gene still overlaps,
            # return the start of the most upstream gene
            if i-j == 0:
                region_start = contig_genes[0].start
                break

            # region starts at first base after  the first 
            # upstream gene that does not overlap
            if contig_genes[i-j].end < contig_genes[i].start:
                region_start = max( contig_genes[i-j].end+1, furthest_gene_end+1 )
                break

        # iterate over downstream genes
        for k in range(1, len(contig_genes)-i+1):

            # if the most downstream gene still overlaps,
            # return the end of the most downstream gene
            if k == len(contig_genes)-i:
                region_end = contig_genes[i+k].end
                break

            # region ends   at last  base before the first 
            # downstream gene that does not overlap
            if contig_genes[i].end < contig_genes[i+k].start:
                region_end = contig_genes[i+k].start-1
                break

    # if any other gene,
    # region starts at first base after  previous gene
    # region ends   at last  base before next     gene
    else:
        region_start = max( contig_genes[i-1].end+1, furthest_gene_end+1 )
        region_end   = contig_genes[i+1].start-1
        
    return (region_start, region_end)


def get_introns_from_bam(bamfile, contig, strand, start, stop, fasta):
    # select good introns supported by reads
    selected_introns = []

    # get contig sequence
    # necessary for intron motif checking
    contig_seq = fasta[contig].seq

    # fetch reads
    # print(f'fetch at {contig=}, {start=}, {stop=}, {strand=}')
    fetched_seqs = bamfile.fetch(contig=contig, start=start, stop=stop)

    # find introns
    ## returned coords are 0-indexed
    ## start = start of intron, end = start of first exon after intron
    ## its thus useful for python slicing
    if strand == '+':
        # introns = bamfile.find_introns( (read for read in fetched_seqs if read.is_forward and read.is_mapped) )
        # introns = bamfile.find_introns( (read for read in fetched_seqs if read.is_read1 and read.is_mapped) )
        introns = bamfile.find_introns( (read for read in fetched_seqs if read.is_mapped and read.has_tag('XS') and read.get_tag('XS') == '+') )
        # allowed_motifs = ["gtag", "gcag", "atac", "atag", "gaag", "acac", "ggag"]
        allowed_motifs = args.splice_sites
    elif strand =='-':
        # introns = bamfile.find_introns( (read for read in fetched_seqs if read.is_reverse and read.is_mapped) )
        # introns = bamfile.find_introns( (read for read in fetched_seqs if read.is_read2 and read.is_mapped) )
        introns = bamfile.find_introns( (read for read in fetched_seqs if read.is_mapped and read.has_tag('XS') and read.get_tag('XS') == '-') )
        # allowed_motifs = ["ctac", "ctgc", "gtat", "ctat", "cttc", "gtgt", "ctcc"]
        # rc_telomere_seq: str = str( Seq(args.telomere_seq).reverse_complement() )
        allowed_motifs = [ str(Seq(site).reverse_complement()) for site in args.splice_sites ]

    # filter away introns that do not need criteria
    for coords, count in introns.items():

        # print(coords, count)
        intron_start, intron_end = coords

        # skip if intron is situated entirely before or after specified region
        # (fetch() also fetches reads that merely overlap with region
        #  some of these reads may indicate introns that are
        #  entirely situated before or after specified region)
        if intron_end < start or stop < intron_start:
            continue

        # skip if intron has less than 5 reads supporting it
        if count < 5:
            continue
        
        # skip if intron is too short or too long
        # intron_len = intron_end - intron_start + 1
        # if not args.min_intron_len < intron_len < args.max_intron_len:
        #     continue

        # skip if intron has non-allowed intron motifs
        donor_site = contig_seq[intron_start:intron_start+2]
        accep_site = contig_seq[intron_end-2:intron_end]
        if f'{donor_site}{accep_site}'.lower() not in allowed_motifs:
            continue

        # skip if intron has poor relative rnaseq support
        if get_intron_rnaseq_support(contig, intron_start, intron_end, bamfile) < args.threshold:
            continue

        # make intron object
        new_intron = gffutils.feature.Feature(
                seqid = contig,
                source = 'fix_genes.py',
                featuretype = 'intron',
                start = intron_start+1,
                end = intron_end,
                score = str(count),
                strand = strand,
                frame = '.',
                attributes = f'ID=new_intron;Parent=placeholder;Name=new_intron'
                )

        # populate list
        # transform intron coords to 1-indexing
        selected_introns.append(new_intron)

    return selected_introns


def select_non_overlapping_introns(introns):
    '''
    Find the set of introns that do not overlap
    Assumes introns are ordered by start coordinate
    '''
    # return nothing if intron list is empty
    if not introns:
        return

    overlapping_introns = []
    # prev_intron_start, prev_intron_end = 0, 100000000
    overlap_start, overlap_end = introns[0].start, introns[0].end
    # iterate over introns one by one
    # assume they are sorted by start, then by end coords
    for i in introns:
        # print(f'{i = }'
        # if current intron overlaps with previous intron,
        # add it to overlapping_introns group
        if overlap_start <= i.start <= overlap_end:
            # print('current intron overlaps with previous intron')
            overlapping_introns.append(i)
            overlap_start = min( overlap_start, i.start)
            overlap_end   = max( overlap_end,   i.end)
        # if current intron does not overlap previous intron,
        # yield the best scored intron from the overlap group
        # and reset the intron group
        elif overlap_end < i.start:
            # print('current intron does not overlap with previous intron')
            best_intron = max( overlapping_introns, key=lambda x : int(x.score) )
            yield best_intron
            overlapping_introns = [i]
            overlap_start, overlap_end = i.start, i.end
        # update start and end coords of previous intron
        # prev_intron_start, prev_intron_end = i_start, i_end
        # overlap_start, overlap_end = max(i_start, overlap_start), i_end
    # for the last intron group, select the best intron
    # after all introns have been iterated over
    if len(introns) > 0:
        best_intron = max( overlapping_introns, key=lambda x : int(x.score) )
        yield best_intron


def splice_seq(seq, introns):
    '''
    Take a sequence string and an iterator or list of introns
    and splice out the introns from the sequence string
    '''

    # do not splice if input intron list is empty
    if len(introns) == 0:
        return seq

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


def get_non_overlapping_orfs(seq: str, strand: str) -> list[tuple[int, int, str, str]]:
    '''
    Given a DNA sequence string, predict the ORFs with orfipy
    and select the biggest ORFs that do not overlap
    '''
    seq = seq.upper()
    # convert strand symbol
    strand = 'f' if strand == '+' else 'r'
    # find orfs in seq
    ## orfs() returns a list of tuples
    all_orfs = orfipy_core.orfs(
        seq,
        minlen=300,
        maxlen=100000,
        starts=['ATG'],
        strand=strand,
        include_stop=True
    )
    non_overlapping_orfs = select_orfs(orfs=all_orfs, overlap='no-overlap')
    return non_overlapping_orfs


def select_orfs(orfs=None, overlap=None) -> list[tuple[int, int, str, str]]:
    selected_orfs = []

    # first sort orfs by length
    orfs.sort(key=lambda x : x[1]-x[0], reverse=True)

    # print(f'all possibly overlapping orfs: {orfs}')
    # then select non-overlapping orfs,
    # preferring the largest ones
    occupied_regions = []
    for i, orf in enumerate(orfs):
        # get coordinates, frame and length of orf
        start, stop = orf[0], orf[1]
        # always select the first orf,
        # i.e. the largest orf
        if i == 0:
            occupied_regions.append([start,stop])
            selected_orfs.append(orf)
            continue
        # check if current orf overlaps with
        # any of the already selected orfs
        for region in occupied_regions:
            if overlap == 'no-overlap':
                # move to next orf if any region overlaps
                if max(region[0], start) < min(region[1], stop):
                    break
            elif overlap == 'partial':
                # move to next orf if this orf is fully within an occupied region
                if region[0] <= start < stop <= region[1]:
                    break

        # this else block is only executed if the above
        # for loop is exited normally (i.e. doesn't break)
        else:
            # if it survived the overlap checks, select orf
            occupied_regions.append([start,stop])
            selected_orfs.append(orf)

    return selected_orfs


# def re_insert_introns(orfs, introns):
def re_insert_introns(orfs: list[tuple[int, int, str, str]],
                      introns
                      ) -> list[list[int, int, str, str]]:
    '''
    Adjust start and end coordinates of predicted ORFs by
    checking which introns overlap with the ORFs
    '''
    # orfs is a list of tuples (0-indexed),
    ## start = 0-indexed,
    ## end   = 0-indexed+1, so actually 1-indexed effectively
    # introns is a list of intron objects (1-indexed)
    # convert to lists of lists so I can edit values
    upd_orfs = [ list(o) for o in orfs ]
    for o in upd_orfs:
        for i in introns:
            # junction is the coord in the spliced seq
            # where the intron is supposed to be
            # it is updated with each intron we visit / re-insert
            junction = i.start - 1 - 1 # -1 to make 0-indexed, -1 to get to junction
            intron_len = i.end - i.start + 1
            # if junction is situated before predicted orf
            if junction < o[0]:
                # print('junction before start')
                # update orf start and end coords
                o[0] += intron_len
                o[1] += intron_len
            # if junction is within the predicted orf
            ## subtract 1 from o[1] to make it 0-indexed,
            ## and thus a 0-index to 0-index comparison
            elif o[0] <= junction < o[1]-1:
                # print('junction within orf')
                # update only the orf end coord
                o[1] += intron_len
            # if junction is after the predicted orf
            elif o[1]-1 < junction:
                # print('junction after orf')
                # move on to the next orf
                break
    return upd_orfs


def create_new_feature(featuretype, orf, templ, gene_number, start=None, end=None, exon_number=None, phase='.'):
    '''
    Create a new feature object.
    Use an existing gene feature as template.
    '''

    # construct gene id, gene name and mrna id from template gene ID and Name
    # gene_id    = re.sub('([0-9]{4})', r'\1-' + str(gene_number), templ['ID'][0])
    gene_id    = re.sub('(_[0-9]+)', r'\1-' + str(gene_number), templ['ID'][0])

    if 'Name' in templ.attributes:
        # gene_name  = re.sub('([0-9]{4})', r'\1-' + str(gene_number), templ['Name'][0])
        gene_name  = re.sub('(_[0-9]+)', r'\1-' + str(gene_number), templ['Name'][0])
    else:
        gene_name = gene_id

    mrna_id    = gene_id.replace('gene','mRNA') + '-T1'

    if featuretype == 'gene':
        start, end, strand, *_  = orf
        attributes = f'ID={gene_id};Name={gene_name}'
    
    elif featuretype == 'mRNA':
        start, end, strand, *_  = orf
        mrna_name  = gene_name.replace('gene','mRNA')
        attributes = f'ID={mrna_id};Parent={gene_id};Name={mrna_name}'

    elif featuretype == 'exon':
        strand     = orf[2]
        exon_id    = f'{mrna_id}.exon{exon_number:02d}'
        attributes = f'ID={exon_id};Parent={mrna_id}'

    elif featuretype == 'CDS':
        strand     = orf[2]
        cds_id     = f'{mrna_id}.CDS'
        attributes = f'ID={cds_id};Parent={mrna_id}'

    new_feature = gffutils.feature.Feature(
            seqid=templ.seqid,
            source='fix_genes.py', 
            featuretype=featuretype, 
            start=start, 
            end=end, 
            score='.',
            strand=strand,  
            frame=str(phase),
            attributes=attributes
    )

    return new_feature


def determine_new_phase(feature, phase):
    possible_phases = [0,2,1]
    codon_remainder = ( len(feature) - phase ) % 3
    phase = possible_phases[codon_remainder]
    return phase



if __name__ == '__main__':
    main(args)
