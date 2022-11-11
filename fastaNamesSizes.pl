#!/usr/bin/perl -w

=head1 SYNOPSIS

fastaNamesSizes.pl - returns a tab list with name and size of sequences 
                     in a fasta file

=head1 USAGE

fastaNamesSizes.pl [-f format] fasta_file1 fasta_file2...

=head2 -f format

Format

=head1 AUTHOR

Lionel Guy (lionel.guy@icm.uu.se)

=cut

use strict;
use Getopt::Std;
use Bio::SeqIO;

our $opt_f;
getopts('f:');

usage() unless @ARGV;

foreach my $file (@ARGV){
    my $count = 0;
    my $sum = 0;
    my $min;
    my $max;
    my $seq_io;
    if ($opt_f){
	$seq_io = Bio::SeqIO->new(-file => "$file", -format => $opt_f);
    }
    else {
	$seq_io = Bio::SeqIO->new(-file => "$file");
    }
    while (my $seq = $seq_io->next_seq){
	my $id = $seq->display_id;
	my $len = $seq->length;
	print "$id\t$len\n";
	# stats
	$count++;
	$sum += $len;
	$min = $len if (!$min || $len < $min);
	$max = $len if (!$max || $len > $max);
    }
    my $average = int($sum/$count+0.5);
    print STDERR "# File $file\n";
    print STDERR "# $count sequences, total length $sum.\n";
    print STDERR "# Minimum len: $min. Max: $max. Average: $average\n";
}
sub usage{
    system("perldoc $0");
    exit;
}
