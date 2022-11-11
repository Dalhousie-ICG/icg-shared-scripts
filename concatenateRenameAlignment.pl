#!/usr/bin/perl -w

=head1 NAME

renameConcatenateAlignment.pl - concatenates gene alignments

=head1 SYNOPSIS

renameConcatenateAlignment.pl [-m|--mapping mapping_file] [-k|--key-values value1,value2] [-f|--format-output format] [-p|--partition-output-file file] [-d|--default-model model] [-g|--gene-list file] [-s|--suffices suffix1,suffix2] fasta_file1 fasta_file2 ...

=head1 DETAILS

Concatenates gene alignments and rename them according to a mapping file. Makes a complete list of all organisms, and prints a complete alignment, filling in with gaps where a sequence is missing. If no mapping file is provided, assumes that ids from sequences represent grouping (organisms).

Optionally, a RAxML-compatible partition file can be produced. In that case, either the same model is applied to all partitions, or a list of partitions (equivalent to genes) is given. 

If the --gene-list option is used, the fasta files must correspond to the genes given in the list, with a .fasta, .fa or .faa suffix, or any other given through --suffixes.

=head1 INPUT

=over

=item B<-m>, B<--mapping> I<mapping_file>

A mapping file that gives mapping for the genes. A tab file.

=item B<-k>, B<--key-values> I<value1,value2>

Which columns of the mapping file give the id and the organism, respectively. Two numbers separated by a comma. By default 3,2.

=item B<-f>, B<--format-output> I<format>

Format for the output file. Fasta by default.

=item B<-p>, B<--partition-output-file> I<file>

Output partition file. This file is compatible with RAxML (see option -q in RAxML). 

=item B<-g>, B<--gene-list> I<file>

List of the genes/partitions to include in the concatenate. Optionally, a second column gives the model to output in the partition file. Tab-separated file.

=item B<-s>, B<--suffices> I<suffix1,suffix2...>

Suffixes for the fasta files, comma-separated. By default: .fasta,.faa,.fa

=item B<-d>, B<--default-model> I<model>

Default model to be applied to all genes/partitions if the --gene-list option is not used or the model lacks in the second column.

=item I<fasta_file1, fasta_file2>

Fasta files, aligned.

=back

=head1 OUTPUT

A concatenated aligment on the standard output.

=head1 AUTHOR

Lionel Guy (lionel.guy@icm.uu.se)

=head1 DATE

Thu Mar  7 15:56:03 CET 2013

=head1 COPYRIGHT AND LICENSE

Copyright 2016 Lionel Guy

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
use File::Basename;
use Bio::SeqIO;
use Bio::Seq;

my %map;
my %seqs;
my %partitions;
my @partitions;


my ($mapping_file, $partition_file, $gene_list_file);
my $suffices = "fasta,faa,faa";
my $keyval = '3,2';
my $output_format = 'fasta';
my $def_model = 'LG';
my $help;

GetOptions(
	   'm|mapping=s' => \$mapping_file,
	   'k|key-values=s' => \$keyval,
	   'f|format-output=s' => \$output_format,
	   'p|partition-output-file=s' => \$partition_file,
	   'd|default-model=s' => \$def_model,
	   'g|gene-list=s' => \$gene_list_file,
	   's|suffices=s' => \$suffices,
	   'h|help' => \$help,
	  );

pod2usage(-exitval => 1, -verbose => 2) if $help;
pod2usage(-exitval => 2, 
	  -message => "No fasta file given\n") unless (@ARGV);

# Suffixes
my @suffixes = split(/,/, $suffices);

# Parse column and organism colums 
my ($id_col, $org_col) = split(/,/, $keyval);
die "No id/org column specified for mapping file\n" 
  unless ($id_col && $org_col && $id_col > 0 && $org_col > 0);

# Read gene list
if ($gene_list_file){
  print STDERR "Reading gene list $gene_list_file\n";
  open GENES, '<', $gene_list_file or die "Cannot open $gene_list_file: $!\n";
  while (<GENES>){
    chomp;
    my ($partition, $model) = split;
    $model = $def_model unless $model;
    $partitions{$partition} = $model;
    push @partitions, $partition;
  }
}

# Read mapping file
if ($mapping_file){
  print STDERR "Reading mapping file $mapping_file\n";
  open MAP, '<', $mapping_file 
    or die "Can't open mapping file $mapping_file\n";
  while (<MAP>){
    my @fs = split;
    my $id = $fs[$id_col-1];
    my $org = $fs[$org_col-1];
    die "No data at colums $id_col or $org_col, line $_\n" 
      unless ($id && $org);
    $map{$id} = $org;
    $seqs{$org}{'n'}++;
  }
  print STDERR "  found " . scalar(keys %map) . " sequences in " . 
    scalar(keys %seqs) . " organisms\n";
}


# Read sequences
print STDERR "Reading " . scalar(@ARGV) . " fasta files\n";
my $nseq = 0;
my $nmissing = 0;
my %seqs_in_fasta;
foreach my $file (@ARGV){
  my ($partition, $path, $suffix) = fileparse($file, @suffixes);
  # Check that the file is in the gene list
  if ($gene_list_file){
    die "File $file have no corresponding partition $partition " . 
    "in $gene_list_file\n" unless ($partitions{$partition});
  } else {
    push @partitions, $partition;
  }
  my $seqin = Bio::SeqIO->new(-file => $file,
			      -format => 'fasta',
			      -idlength => 30);
  my $length;
  while (my $seq = $seqin->next_seq){
    $nseq++;
    $length = $seq->length unless $length;
    die "Sequence " . $seq->id . "has a different length (" . $seq->length .
      ") as the rest ($length)\n" unless ($length == $seq->length);
    my $org;
    if ($mapping_file){
      $org = $map{$seq->id};
      die "No mapping for " . $seq->id . " in file $file\n" unless $org;
    } else {
      $org = $seq->id;
      $seqs{$org}{'n'}++;
    }
    $seqs_in_fasta{$partition}{$org}{'seq'} = $seq;
    }
  $seqs_in_fasta{$partition}{'length'} = $length;
}

# Go through partitions, build concatenates
my $position = 1;
if ($partition_file){
  open PARTITION, '>', $partition_file or die 
    "Cannot open $partition_file for writing\n";
}
foreach my $partition (@partitions){
  # Check that we have data for all partitions
  die "No fasta file for partition $partition\n" 
    unless $seqs_in_fasta{$partition};
  # Loop through organisms
  foreach my $org (keys %seqs){
    if ($seqs_in_fasta{$partition}{$org}){
      $seqs{$org}{'seq'} .= $seqs_in_fasta{$partition}{$org}{'seq'}->seq;
    } else {
      $seqs{$org}{'seq'} .= "-" x $seqs_in_fasta{$partition}{'length'};
      $nmissing++;
    }
  }
  # Print partition file if requested
  if ($partition_file){
    my $model = $def_model;
    $model = $partitions{$partition} if ($partitions{$partition});
    print PARTITION "$model, $partition = $position" . "-" . 
      ($position + $seqs_in_fasta{$partition}{'length'} - 1) . "\n";
  }
  $position += $seqs_in_fasta{$partition}{'length'};
}
# Report on missing data
print STDERR sprintf "  %d alignments in %d species\n  %.1f missing\n", 
    scalar(@ARGV), scalar(keys %seqs), $nmissing/$nseq*100;

# Print fasta
print STDERR "Printing concatenated alignment\n";
my $seqout = Bio::SeqIO->new(-fh => \*STDOUT,
			     -format => $output_format);
foreach my $org (sort keys %seqs){
  my $newseq = Bio::Seq->new(-seq => $seqs{$org}{'seq'},
			     -id => $org);
  $seqout->write_seq($newseq);
}

