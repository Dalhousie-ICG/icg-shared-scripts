#!/usr/bin/perl -w
use strict;
use Bio::DB::EUtilities;
use Bio::SeqIO;

my ($list,$out) = @ARGV;
my @accessions;

open IN, "<$list";
while (<IN>) {
    chomp;
    push @accessions, $_;
}
close IN;

# print join "\n", @accessions;

my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
                                       -db      => 'nucleotide',
                                       -rettype => 'gbwithparts',
                                       -email   => 'mymail@foo.bar',
                                       -id      => \@accessions);
 
my $file = "$out";
$factory->get_Response(-file => $file);
