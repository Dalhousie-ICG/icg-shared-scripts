#!/usr/bin/env python

import gffutils
import pysam
from Bio import SeqIO
import orfipy_core
import argparse
from statistics import mean
import re

# make command line interface
parser = argparse.ArgumentParser(
        description=
        """
        Checks whether predicted introns are supported by RNAseq data,
        and removes them if they are not.

        Strips away all predicted genes that have false introns, and
        replaces them with a set of non-overlapping ORFs. If a
        predicted gene also has one or more true introns, they will
        be re-inserted into the new ORFs.

        If a region has no read mapping at all, and an intron was
        predicted there, the script will keep the intron.
        """)
parser.add_argument(
    "-g", "--gff3",
    dest="gff3_file",
    metavar="GFF3_FILE",
    required=True,
    help="Input GFF3 file")
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
args = parser.parse_args()

def add_rnaseq_supports(introns, bamfile):
    '''
    Calculate rnaseq supports for introns and add 
    as new attribute to intron objects.
    RNAseq support = n reads that span intron / all reads
    '''
    introns_w_support = []
    for i in introns:
        # start list of spliced / all aligned reads ratio's
        ratios = []
        # i.start and i.end are 1-indexed
        # pysam works with 0-indexing
        # so to get coverage for position x (1-indexed),
        # we need to ask for position x-1 (0-indexed)
        for col in bamfile.pileup(contig=i.seqid, start=i.start-1, stop=i.end, truncate=True):
            num_aligned = col.get_num_aligned()
            # skip col if no reads are aligned
            if num_aligned == 0:
                continue
            # start spliced alignment counter
            spliced = 0
            # iterate over aligned reads at col
            for read in col.pileups:
                # add 1 to spliced alignment counter if read has 'N' here
                # in the CIGAR string. This indicates the read is spliced here
                if read.is_refskip:
                    spliced += 1
            # calculate ratio (spliced / all aligned reads)
            ratio = spliced / num_aligned
            # append ratio to list of ratios
            ratios.append(ratio)
        # if no reads map to this intron at all,
        # assign ratio 'NA'
        if len(ratios) > 0:
            i.attributes['ratio'] = f'{mean(ratios):.2f}'
        else:
            i.attributes['ratio'] = '1.00'
        introns_w_support.append(i)
    return introns_w_support

def get_region(db, gene, fasta):
    '''
    Get the start and end coordinates of a sequence region,
    as well as its DNA sequence. The region stretches from 
    the first base after the previous gene until
    the first base before the next gene
    '''
    # deconstruct ID attribute of gene
    gene_id = gene['ID'][0]
    gene_number = int( gene_id.split('gene')[1] )
    gene_contig = str( gene_id.split('gene')[0] )
    # calculate total number of genes on contig
    genes_on_contig = list( db.region(seqid=gene.seqid,start=1,featuretype='gene') )

    # if first gene on contig,
    # region starts at first base of     first gene
    # region ends   at last  base before next  gene
    if gene_number == 1:
        next_gene_id = gene_contig + "gene%04d" % (gene_number + 1)
        region_start = db[gene_id].start
        region_end   = db[next_gene_id].start - 1
    # if last gene on contig,
    # region starts at first base after previous gene
    # region ends   at last  base of    current  gene
    elif gene_number == len(genes_on_contig):
        prev_gene_id = gene_contig + "gene%04d" % (gene_number - 1)
        region_start = db[prev_gene_id].end + 1
        region_end   = db[gene_id].end
    # if any other gene,
    # region starts at first base after  previous gene
    # region ends   at last  base before next     gene
    else:
        prev_gene_id = gene_contig + "gene%04d" % (gene_number - 1)
        next_gene_id = gene_contig + "gene%04d" % (gene_number + 1)
        region_start = db[prev_gene_id].end   + 1
        region_end   = db[next_gene_id].start - 1
    
    # retrieve contig DNA sequence
    contig_seq = fasta[gene.seqid].seq
    # slice out region sequence
    ## python slicing is 0-indexed, so substract 1
    ## from 1-indexed region_start value
    region_seq = str( contig_seq[region_start-1:region_end] )
    return ( region_start, region_end, region_seq )

def select_non_overlapping_introns(introns):
    '''
    Find the set of introns that do not overlap
    '''
    selection = []
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
            best_intron = max( overlapping_introns, key=lambda x : float(x['ratio'][0]) )
            selection.append(best_intron)
            overlapping_introns = [i]
        # update start and end coords of previous intron
        prev_intron_start, prev_intron_end = i.start, i.end
    # for the last intron group, select the best intron
    # after all introns have been iterated over
    if len(introns) > 0:
        best_intron = max( overlapping_introns, key=lambda x : float(x['ratio'][0]) )
        selection.append(best_intron)
    return selection

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

def get_non_overlapping_orfs(seq):
    '''
    Given a DNA sequence string, predict the ORFs with orfipy
    and select the biggest ORFs that do not overlap
    '''
    non_overlapping_orfs = []
    seq = seq.upper()
    # # convert strand symbol
    # if strand == '+':
    #     strand = 'f'
    # elif strand == '-':
    #     strand = 'r'
    # find orfs in seq
    ## orfs() returns a list of tuples
    all_orfs = orfipy_core.orfs(
        seq,
        minlen=300,
        maxlen=10000,
        starts=['ATG'],
        strand='b',
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
            non_overlapping_orfs.append(orf)
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
            non_overlapping_orfs.append(orf)
    return non_overlapping_orfs

def adjust_orfs_w_introns(orfs, introns, region_start):
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
            elif o[0] < junction < o[1]:
                # update only the orf end coord
                o[1] += intron_len
            # if junction is after the predicted orf
            elif o[1] < junction:
                # move on to the next orf
                break
        # transform back to genome space
        o[0] = o[0] + region_start
        o[1] = o[1] + region_start - 1
    # sort updated orfs by start position
    upd_orfs.sort(key=lambda x : x[0])
    return upd_orfs

def create_new_feature(featuretype, orf, templ, gene_number, start=None, end=None, exon_number=None, phase='.'):
    '''
    Create a new feature object.
    Use an existing gene feature as template.
    '''

    # construct gene id, gene name and mrna id from template gene ID and Name
    gene_id    = re.sub('([0-9]{4})', r'\1-' + str(gene_number), templ['ID'][0])
    gene_name  = re.sub('([0-9]{4})', r'\1-' + str(gene_number), templ['Name'][0])
    mrna_id    = gene_id.replace('gene','mRNA')

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
            source='script', 
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


# load gff file with introns explicitly defined
db = gffutils.create_db(args.gff3_file, ':memory:', merge_strategy='create_unique')

# load bam file
bamfile = pysam.AlignmentFile(args.bam_file, mode='rb')

# load fasta
fasta = SeqIO.index(args.fasta_file, 'fasta')

# main code
false_genes = []
new_features = []
for g in db.features_of_type('gene'):
    # print()
    # print(g.id)

    # get introns belonging to the gene
    introns = db.children(g, featuretype='intron')
    # add rnaseq support ratios to intron objects
    introns_w_support = add_rnaseq_supports(introns, bamfile)

    # sort introns to true and false introns
    true_introns  = [ i for i in introns_w_support if float(i['ratio'][0]) >= args.threshold ]
    false_introns = [ i for i in introns_w_support if float(i['ratio'][0]) <  args.threshold ]
    # delete false introns from database
    db.delete(false_introns)

    # if any introns do not have sufficient support,
    # re-predict genes in this area
    if len(false_introns) > 0:
        # print()
        # print(g)

        # mark gene for deletion
        false_genes.append(g)

        # region = previous gene end + 1 and next gene start - 1
        # get region coords and region sequence
        region_start, region_end, region_seq = get_region(db, g, fasta)

        # ensure introns do not overlap
        true_introns = select_non_overlapping_introns(true_introns)

        # transform intron coords to region_seq space 
        for i in true_introns:
            i.start = i.start-region_start+1 
            i.end   = i.end  -region_start+1 

        # splice introns from region_seq
        spliced_region_seq  = splice_seq(region_seq, true_introns)
        # print('spliced_seq:', spliced_region_seq)

        # predict orfs
        new_orfs = get_non_overlapping_orfs(spliced_region_seq)

        # skip to next gene if no orfs were found
        if len(new_orfs) == 0: continue

        # adjust orf start and stop by checking for overlapping introns
        # also moves orfs back to genome space
        adj_orfs = adjust_orfs_w_introns(new_orfs, true_introns, region_start)
        # print(adj_orfs)
        
        # transform intron coords back to genome space 
        for i in true_introns:
            i.start = i.start+region_start-1 
            i.end   = i.end  +region_start-1 

        # loop over orfs and create new gene and mRNA for each new orf
        n = 1 # gene counter
        for orf in adj_orfs:
            # print()
            # print(orf)
            # create new gene
            new_gene = create_new_feature('gene', orf, templ=g, gene_number=n)
            # print(new_gene)

            # create new mRNA
            new_mrna = create_new_feature('mRNA', orf, templ=g, gene_number=n)
            # print(new_mrna)

            # add to new_features
            new_features.extend( [new_gene,new_mrna] )

            # loop over introns and create new exons, CDSs and introns along the way
            phase = 0
            m = 1 # exon counter

            # if orf is on + strand
            if orf[2] == '+':
                moving_start = orf[0]
                for i in true_introns:
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
                        new_exon = create_new_feature('exon', orf, templ=g, gene_number=n, start=exon_start, end=exon_end, exon_number=m)
                        # print(new_exon)
                        # create new CDS
                        new_cds = create_new_feature('CDS', orf, templ=g, gene_number=n, start=exon_start, end=exon_end, phase=phase)

                        # add to new_features
                        new_features.extend( [new_exon,new_cds] )

                        # print(new_cds)
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
                new_exon = create_new_feature('exon', orf, templ=g, gene_number=n, start=exon_start, end=exon_end, exon_number=m)
                # print(new_exon)
                new_cds = create_new_feature('CDS', orf, templ=g, gene_number=n, start=exon_start, end=exon_end, phase=phase)
                # print(new_cds)

                # add to new_features
                new_features.extend( [new_exon,new_cds] )

            elif orf[2] == '-':
                true_introns.reverse()
                moving_end = orf[1]
                for i in true_introns:
                    # if intron after gene
                    if orf[1] < i.start:
                        continue
                    # if intron within gene
                    elif orf[0] < i.end < orf[1]:
                        # create new exon
                        exon_start  = i.end + 1
                        exon_end    = moving_end
                        new_exon = create_new_feature('exon', orf, templ=g, gene_number=n, start=exon_start, end=exon_end, exon_number=m)
                        # print(new_exon)
                        # create new CDS
                        new_cds = create_new_feature('CDS', orf, templ=g, gene_number=n, start=exon_start, end=exon_end, phase=phase)
                        # print(new_cds)

                        # add to new_features
                        new_features.extend( [new_exon,new_cds] )

                        # update phase, moving end and counter for next exon
                        phase = determine_new_phase(new_exon, phase)
                        moving_end = i.start - 1
                        m += 1
                    # if intron before gene
                    elif i.end < orf[0]:
                        true_introns.reverse()
                        break
                # the last exon and CDS
                # or first exon and CDS if there were no introns
                exon_start  = orf[0] 
                exon_end    = moving_end
                new_exon = create_new_feature('exon', orf, templ=g, gene_number=n, start=exon_start, end=exon_end, exon_number=m)
                # print(new_exon)
                new_cds = create_new_feature('CDS', orf, templ=g, gene_number=n, start=exon_start, end=exon_end, phase=phase)
                # print(new_cds)

                # add to new_features
                new_features.extend( [new_exon,new_cds] )

            # update gene counter
            n += 1


# delete the old gene with false introns
# and delete all its children
for g in false_genes:
    for f in db.children(g, featuretype=('exon','intron','CDS','mRNA')):
        db.delete(f)
    db.delete(g)

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
for f in db.all_features(order_by=('seqid','start','featuretype'), reverse=True):
    if f.featuretype == 'mRNA':
        print()
        print(f)
    else:
        print(f)
