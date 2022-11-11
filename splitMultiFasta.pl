#!/usr/bin/perl

=head1 NAME

    splitMultiFasta.pl - splits a single multifasta file into multiple singlefasta files

=head1 USAGE

    splitMultiFasta.pl -i|--multifasta <multifasta> -o|--outdir <outdir>

=head1 DESCRIPTION

    Splits a single multifasta file into multiple singlefasta files.
    Concatenates fasta id and fasta description with _
    Also converts any / into _ in id or description
    
=head1 AUTHOR

    Joran Martijn (joran.martijn@icm.uu.se)
    
=cut

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Pod::Usage;

my ($multifasta,$outdir);
GetOptions("i|multifasta=s" => \$multifasta,
	   "o|outdir=s"     => \$outdir,
);
pod2usage (-msg => "Not enough arguments") unless $outdir && $multifasta;

my $in  = Bio::SeqIO->new(-format => 'fasta',
			  -file   => $multifasta);

while (my $seq = $in->next_seq) {

    # retrieve seqname
    my $id   = $seq->id;
    my $desc = $seq->desc;
    
    # make filename
    my $singlefasta = $id . '_' . $desc . '.fasta';
    $singlefasta =~ s/[\s+\/]/_/g;
    print $singlefasta, "\n";

    # make outfile
    my $out = Bio::SeqIO->new(-format => 'fasta',
    			      -file   => ">$outdir/$singlefasta");
    $out->write_seq($seq);
}
