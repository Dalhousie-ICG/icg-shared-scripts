#!/usr/bin/perl -w

=head1 NAME

    getCodingDensity.pl - estimates coding density of a sequence using its genbank file

=head1 USAGE

    getCodingDensity.pl -g|gbk [.gbk] -f|features [comma separated features] > [report.out]

=head1 DESCRIPTION

    Estimates the coding density of a sequence using its genbank file.
    Works also when there are multiple contigs within the genbank file.
    Incorporates also ori spanning features.

=head1 EXAMPLE

    getCodingDensity.pl -g Holospora_undulata_HU1.gbk -f CDS,rRNA,tRNA > Holospora_codingDensity.report

=head1 TODO

    Script is still quite ineffecient.
    Will load regions that code for multiple coding features (for example overlapping CDS) multiple times in the hash.
    This should be in principle unnecessary.

=head1 AUTHOR

    Joran Martijn (joran.martijn@icm.uu.se)

=cut

use strict;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

pod2usage (-msg => "Not enough arguments") unless @ARGV == 4;

my ($gbk,$fts);
GetOptions('g|gbk=s'      => \$gbk,
	   'f|features=s' => \$fts);

my $in = Bio::SeqIO->new(-file   => $gbk,
			 -format => 'genbank');

my @fts = split ',', $fts;

my $total_coding = 0;
my $total_genome_size = 0;
my $contig_count = 0;

while (my $seq = $in->next_seq) {

    # declare data structures
    my %coding_sites;
    my $contig_size;

    # add 1 to contig count
    $contig_count++;

    # loop over features of first contig
    for my $ft ($seq->get_SeqFeatures) {

	# get contig size
	$contig_size = $ft->end() if ($ft->primary_tag eq 'source');

	# skip all non-coding features
	next unless ($ft->primary_tag ~~ @fts);

	## count coding sites
	# if CDS spans ori
	if ( $ft->location->isa('Bio::Location::SplitLocationI') ) {
	    for my $loc ( $ft->location->sub_Location ) {
		# print $loc->start . ".." . $loc->end . "\n";
		my ($start, $end) = ($loc->start, $loc->end);
		for (my $i = $start; $i <= $end; $i++) {
		    $coding_sites{$i}++;
		}

	    }
	}

	# all other CDSs
	else {
	    # print $start, "\t", $end, "\n";
	    my ($start, $end) = ($ft->start, $ft->end);
	    for (my $i = $start; $i <= $end; $i++) {
		$coding_sites{$i}++;
	    }
	}

    }

    # add number of coding sites of first contig to total count
    my $coding = scalar keys %coding_sites;
    $total_coding += $coding;

    # add contig size to total assembly/genome size
    $total_genome_size += $contig_size;

}

print "NUMBER OF CONTIGS: ", "\t", $contig_count, "\n";
print "CODING: ", "\t", $total_coding, "\n";
print "GENOME SIZE: ", "\t", $total_genome_size, "\n";
print "CODING DENSITY: ", "\t"; printf("%.3f", $total_coding / $total_genome_size); print "\n";
