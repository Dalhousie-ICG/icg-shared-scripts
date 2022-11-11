#!/usr/bin/perl -w

=head1 NAME

    getIntergenicSpace.pl - Estimates average intergenic space length of a given genome using its genbank file.

=head1 USAGE

    getIntergenicSpace.pl -g|gbk [.gbk] -f|features [comma separated features] > [report.out]

=head1 DESCRIPTION

    Estimates the average intergenic space length of a given genome.
    Works also if genome is in draft state (that is: multiple contigs).
    Does not incorporate ORI spanning region

    In case of multiple contigs, the contig edges that are not coding will not be taken into account for calculating the intergenic regions.
    So only intergenic regions are used of which we know how big they are.
    Non coding contig edges we are not sure how big they are and can thus not be used in the calculation.

    In case of overlapping coding regions the size of the overlap will NOT be substracted from the total intergenic region length.

=head1 EXAMPLE

    getIntergenicSpace.pl -g Holospora_undulata_HU1.gbk -f CDS,rRNA,tRNA > Holospora_intergenicSpace.report

=head1 AUTHOR

    Joran Martijn (joran.martijn@icm.uu.se) with help of Guilleaume Reboul.

=cut

use strict;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

pod2usage (-msg => "Not enough arguments")  unless @ARGV >= 4;

my ($gbk,$fts);
GetOptions('g|gbk=s'      => \$gbk,
	   'f|features=s' => \$fts);

my $in = Bio::SeqIO->new(-file   => $gbk,
			 -format => 'genbank');

my @fts = split ',', $fts;

my $total_coding = 0;
my $trim_genome_size = 0;
my $total_genome_size = 0;
my $total_ft_count = 0;
my $total_intergenic_region_count = 0;

# loop over contigs
while ( my $seq = $in->next_seq) {

    my %coding_sites;
    my $ft_count = 0;
    my ($init_start, $init_end);
    my $trim_ctg_size;
    my $contig_size;

    # skip contigs without features (only feature is 'source')
    next if ($seq->feature_count == 1);

    # loop over features
    for my $ft ($seq->get_SeqFeatures) {

	# get contig size
	$contig_size = $ft->end() if ($ft->primary_tag eq 'source');

	# skip all non-coding features
	next unless ($ft->primary_tag ~~ @fts);

	# count features
	$ft_count++;
	$total_ft_count++;

	# count coding sites
	my ($coding_start, $coding_end) = ($ft->start, $ft->end);
	for (my $i = $coding_start; $i <= $coding_end; $i++) {
		$coding_sites{$i}++;
	}

	# store coding feature start and end
	$init_start = $ft->start if ($ft_count == 1);
	$init_end   = $ft->end;

#	print $ft_count, "\t", $init_start, "\t", $init_end, "\n";
    }

    # Add number of coding sites to total count
    $total_coding += scalar keys %coding_sites;

    # Add contig size to total genome size
     $total_genome_size += $contig_size;

    # Get the coordinates of the first coding region start
    # and the coordinates of the last  coding region end
    my $start = $init_start;
    my $end   = $init_end;

    # Calculate "trimmed" contig length
    $trim_ctg_size = $end - $start;

    # Add "trimmed" contig size to total "trimmed" genome size
    $trim_genome_size += $trim_ctg_size;

    # Add intergenic region count of contig to total intergenic region count
    my $intergenic_region_count = $ft_count - 1;
    $total_intergenic_region_count += $intergenic_region_count;
}

# Calculate total size of non-coding regions in between coding regions
my $total_non_coding = $trim_genome_size - $total_coding;

# Calculate average intergenic region size
my $avg_intergenic_region_size = $total_non_coding / $total_intergenic_region_count;

# Print results
print "TOTAL CODING:", "\t", $total_coding, "\n";
print "TOTAL NONCODING:", "\t", $total_non_coding, "\n";
print "TOTAL TRIMMED GENOME LENGTH:", "\t", $trim_genome_size, "\n";
print "TOTAL GENOME LENGTH:", "\t", $total_genome_size,"\n";
print "AVERAGE INTERGENIC SPACE LENGTH:", "\t"; printf("%.3f", $avg_intergenic_region_size); print "\n";
