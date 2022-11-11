#!/usr/bin/perl -w

=head1 SYNOPSIS

parseSplitCounts.pl - Count predefined splits in a number of split count files

=head1 USAGE

parseSplitCounts.pl -t|--test-split-file splits [-c|--hss-cutoff int] [-s|--suffix str] split_count_files ...

=head1 INPUT

=head2 -t|--test-split-file splits

A file describing the splits to test, one per row, with following format: name, a tab, the split (each species separated by a comma) separated by a space.
E.g. 
a     x,y,z w,v,u
b     x,z y,w,v,u

=head2 split_count_files ...

Split count files, with the number of times the split is supported, a space, both sides of the splits (see above). E.g.:

92    x,y,z w,v,u
34    x,z y,w,v,u

=head2 [-c|--hss-cutoff int]

If set, returns, at the end of the report, the number of highly-supported splits (HSSs), splits present this number or more times.

=head2 [-s|--suffix str]

Suffix to remove from splitcount files, to parse levels. Otherwise, levels are taken from file names. '.splitcount' by default.

=head1 OUTPUT

A tab file with the tested splits in a row and the levels of the splits (...) as columns.

=head1 AUTHOR

Lionel Guy (lionel.guy@icm.uu.se)

=head1 DATE

<date>

=cut

# Libraries
use strict;
use Getopt::Long;
use Sort::Naturally;
use File::Basename;

# Options 
my $test_split_file;
my $suffix = '.splitcount';
my $hssc;

GetOptions(
	   't|test-split-file=s' => \$test_split_file,
	   's|suffix=s'          => \$suffix,
	   'c|hss-cutoff=s'      => \$hssc,
);
usage() unless (@ARGV && $test_split_file);

my %tested_splits;
my @tested_splits_names;
my %split_count;
my %hss_count;
my %levels;

# Parse file that contains tested splits
open TESTED, '<', $test_split_file or die "Could not open $test_split_file\n";
while (<TESTED>){
    chomp;
    next if /^#/;
    my @a = split;
    $tested_splits{$a[1]} = $a[0];
    $tested_splits{$a[2]} = $a[0];
    push @tested_splits_names, $a[0];
}

# foreach (sort keys %tested_splits){
#     print $tested_splits{$_}, "-> $_\n";
# }
# die;

# Parse split count files
foreach (@ARGV){
  chomp;
  # Parse level. Not very good way, but anyway...
  my $lvl = $_;
  $lvl = basename($lvl, ($suffix)) if $suffix;
  #if (/(\d+)/) {$lvl = $1} else {die "No digits found in $_\n";}
  $levels{$lvl}++;
  # Read split count file
  my ($split_count_ref, $nhss) = parseCountFile($_, \%tested_splits, $hssc);
  $hss_count{$lvl} = $nhss;
  foreach (@tested_splits_names){
    if ($split_count_ref->{$_}){
      $split_count{$lvl}{$_} = $split_count_ref->{$_}
    }
    else {
      $split_count{$lvl}{$_} = 0;
    }
  }
}

# Print results: rows are tested splits, columns are levels
# Print title;
print "split";
foreach my $lvl (nsort (keys %levels)){
  print "\t$lvl";
}
print "\n";
# Print the rest
foreach my $split (@tested_splits_names){
  print "$split";
  foreach my $lvl (nsort (keys %levels)){
    if ($split_count{$lvl}{$split}){
      print "\t" . $split_count{$lvl}{$split};
    }
    else {
      print "\t0";
    }
  }
  print "\n";
}
# Print number hss
if ($hssc){
  print "HSS_count";
  foreach my $lvl (nsort (keys %levels)){
    print "\t", $hss_count{$lvl};
  }
}



# Parse a single count file. Returns
sub parseCountFile {
  my ($file, $tested_splits_ref, $hssc) = @_;
  my %seen;
  my $nhss = 0;
  open SPLITS, '<', $file or die;
  while (<SPLITS>){
    chomp;
    my ($n, $s1, $s2) = split;
    # Count highly-supported splits if required
    $nhss++ if ($hssc && $n >= $hssc);
    # Test if any of the splits is present in the splits to test
    if ($tested_splits_ref->{$s1}){
      die "Found only one side of the split in file $file, split $_\n" 
	unless $tested_splits_ref->{$s2};
      die "Non-matching split sides($s1 <--> $s2) in file $file: $_\n"
	unless ($tested_splits_ref->{$s1} eq $tested_splits_ref->{$s2});
      $seen{$tested_splits_ref->{$s1}} += $n;
    }
  }
  return \%seen, $nhss;
}

sub usage{
  system("perldoc $0");
  exit;
}
