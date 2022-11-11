#!/usr/bin/perl -w

=head1 SYNOPSIS

aa_recoding.pl - Recodes amino acid sequences, for example in Dayhoff categories

=head1 USAGE

aa_recoding.pl -a|--alphabet <string> [-r|--recoding <string>] <fasta_file>

=head1 INPUT

=head2 -a|--alphabet <string>

Either a string giving the letters to replace (e.g. ARNDCQEGHILKMFPSTWYVX), or one of "dayhoff6", "dayhoff4", "sr4" or "hp". In the former case, the recoding argument is mandatory. In the latter case, the recoding is done according to the following scheme:

=over

=item dayhoff6

A,G,P,S,T = A;  D,E,N,Q = D; H,K,R = H; F,Y,W = F; I,L,M,V = I; C = C

=item dayhoff4

A,G,P,S,T = A; D,E,N,Q = D; H,K,R = H; F,Y,W,I,L,M,V = F; C = ?

=item hp

A,C,F,G,I,L,M,V,W = A; D,E,H,K,N,P,Q,R,S,T,Y = D

=item sr4

A,G,N,P,S,T = A; C,H,W,Y = C; D,E,K,Q,R = G; F,I,L,M,V = T

=back 

=head2 [-r|--recoding <string>]

A string of the same length as --alphabet. Letter in position 1 in the alphabet will be recoded into the letter in position 1 in this argument, etc. 

=head1 OUTPUT

A fasta file.

=head1 AUTHOR

Lionel Guy (lionel.guy@icm.uu.se)
Joran Martijn (joran.martijn@icm.uu.se) [added SR4 Recoding]

=head1 DATE

Fri May 16 09:57:55 CEST 2014

=head1 COPYRIGHT

Copyright (c) 2014 Lionel Guy

=head1 LICENSE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

# libraries
use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO;

my %alphabets = (
    dayhoff6 => {
	alphabet => 'AGPSTDENQHKRFYWILMVC',
	recoding => 'AAAAADDDDHHHFFFIIIIC',
    },
    dayhoff4 => {
	alphabet => 'AGPSTDENQHKRFYWILMVC',
	recoding => 'AAAAADDDDHHHFFFFFFF?',
    },
    hp => {
	alphabet => 'ACFGILMVWDEHKNPQRSTY',
	recoding => 'AAAAAAAAADDDDDDDDDDD',
    },
    sr4 => {
        alphabet => 'AGNPSTCHWYDEKQRFILMV',
        recoding => 'AAAAAACCCCGGGGGTTTTT',
    },
);


my $alphabet;
my $recoding;
my $help;

GetOptions(
    'a|alphabet=s' => \$alphabet,
    'r|recoding=s'  => \$recoding,
    'h|help' => \$help,
);

pod2usage(-exitval => 1, -verbose => 2) if $help;
pod2usage(-exitval => 2, -message => "No fasta input\n") unless (@ARGV);
pod2usage(-exitval => 2, -verbose => 2,
	  -message => "No alphabet input\n") unless ($alphabet);

if ($alphabets{$alphabet}){
    $recoding = $alphabets{$alphabet}{'recoding'};
    $alphabet = $alphabets{$alphabet}{'alphabet'};
}
pod2usage(-exitval => 1, 
	  -message => "Alphabet $alphabet and recoding $recoding have " . 
	      "different lengths")
    unless (length($alphabet) == length($recoding));


print STDERR "Recoding $alphabet\n";
print STDERR "    with $recoding\n";

my $fasta_in = Bio::SeqIO->new(-file => $ARGV[0],
			       -format => 'fasta');
my $fasta_out = Bio::SeqIO->new(-fh => \*STDOUT,
				-format => 'fasta');
while (my $seq = $fasta_in->next_seq){
    my $s = $seq->seq;
    #print "$s\n";
    $_ = $s;
    eval "tr/$alphabet/$recoding/, 1" or die $@;
    $s = $_;
    #print "$s\n";
    $seq->seq($s);
    $fasta_out->write_seq($seq);
}


