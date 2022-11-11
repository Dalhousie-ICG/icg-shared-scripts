#!/usr/bin/env perl
package alignment_pruner;
=head1 NAME

alingment_pruner.pl

=head1 SYNOPSIS

    alignment_pruner.pl --file alignment.fasta --gap_threshold 10 > pruned.fasta

=head1 DESCRIPTION

alignment_pruner.pl removes unconserved or gappy columns from an alignment
according to criteria specified by the user.

this version contains an extra line of code (line 1035) compared to the
previous published version (Martijn et al 2018) to show the delta chi-square 
value per alignment site when invoking --chi2_prune Plot

The chi2 and overview functions currently only work for amino-acids.

=head1 OPTIONS

=cut

use Moose;
use MooseX::Types::Path::Class;
with 'MooseX::Getopt';

use feature ':5.10';
#no warnings 'experimental';
use experimental 'smartmatch';

use List::MoreUtils qw( uniq minmax );
use List::Util qw( max sum );
use Bio::AlignIO;
use Bio::Matrix::IO;
use Bio::Tools::CodonTable;
use GD;
use Statistics::Descriptive;

#use Benchmark qw/ :all :hireswallclock /;
#use Data::Dump qw( dd );

=head2 --file <file>

The file that contains the alignment.

=cut

has 'file' => (
    is => 'ro',
    isa => 'Path::Class::File',
    required => 1,
    coerce => 1,
    documentation => 'File containing alignment',
);

=head2 --format <format>

The alignment format, default is fasta. Can use any format supported by
bioperl.

=cut

has 'format' => (
    is  => 'ro',
    isa => 'Str',
    required => '',
    default => 'fasta',
    documentation => 'Alignment format, default is fasta',
);

=head2 --subset <subset>

Specify a subset of the sequences that all calculations should be based on.
This is usefull for only taking the ingroup into account for example. The
argument has to be a regular expression that matches the names of the sequences
to include.

    --subset '^INGROUP'

    --subset '^(Snuffe|Snuffa|Peter)$'

=cut

has 'subset' => (
    is => 'ro',
    isa => 'Str',
    required => '',
    default => '',
    documentation => 'The subset of sequences to be used for the analysis, argument is a regexp that should match the id of the sequences',
);

=head2 --chi2_test

This will run a chi2 test on you alignment and show the result as a table.
The table contains 4 columns:

    1. The number of the sequence in the alignment
    2. The name of the sequence
    3. The chi2 statistic
    4. A star if it is significantly deviating

Note that if you have specified a subset, the expected distribution of
aminoacids is based on this set, but the script will calculate the chi2
test on all sequences.

=cut

has 'chi2_test' => (
    is => 'ro',
    isa => 'Bool',
    required => '',
    default => '',
    documentation => 'Show chi2 statistics for the taxa of the alignment',
);

=head2 --bowker_symmetry_test

This will run bowkers symmetry test on you alignment and show the result as a
table.

=cut

has 'bowker_symmetry_test' => (
    is => 'ro',
    isa => 'Bool',
    required => '',
    default => '',
    documentation => 'Run the Bowker symmetry test',
);

=head2 --aminogc

Calculate the aminoGC content of the different taxa.

=cut

has 'aminogc' => (
    is => 'ro',
    isa => 'Bool',
    required => '',
    default => '',
    documentation => 'Calculate aminogc',
);

=head2 --chi2_remove_taxa

This will run a chi2 test on you alignment and remove all taxa that are
significantly deviating.

=cut

has 'chi2_remove_taxa' => (
    is => 'ro',
    isa => 'Bool',
    required => '',
    default => '',
    documentation => 'Remove taxa that fail the chi2 test',
);

=head2 --remove_columns <columns>

Specify a list of columns to remove from the alignment. The first column is
column 0. It should be a comma separated string of ranges. Example:

    0-5,22-30

=cut

has 'remove_columns' => (
    is => 'ro',
    isa => 'Str',
    required => '',
    default => '',
    documentation => 'Remove columns, specified by a commaseparated string (4-10,88,100-140)',
);

=head2 --gap_threshold <threshold>

The threshold used for removing gapped positions, either the maximal number of
sequences that are allowed to have gaps or the fraction of sequences allowed to
have gaps, either as percent or as a number between 0 and 1. Examples:

    --gap_threshold 10   # remove columns with gaps in more than 10 sequences.

    --gap_threshold 10%  # remove columns with gaps in more than 10% of the sequences.

    --gap_threshold .1   # remove columns with gaps in more than 10% of the sequences.

=cut

has 'gap_threshold' => (
    is => 'ro',
    isa => 'Str',
    required => '',
    default => '',
    documentation => 'The threshold used for removing gapped positions, either a number specifying the minimal number of sequences or %'
);

=head2 --conserved_threshold <threshold>

The threshold used for removing unconserved positions. Conservation is
calculated as the number of times the most frequent aminoacid appears in an
alignment column. Specify in the same way as for the --gap_threshold. This
option is probably not that useful. Example:

    --conserved_threshold 10

    --conserved_threshold 10%

    --conserved_threshold .1

=cut

has 'conserved_threshold' => (
    is => 'ro',
    isa => 'Str',
    required => '',
    default => '',
    documentation => 'The threshold used for removing unconserved positions, same as gap_threshold',
);

=head2 --chi2_prune <half|n#|f#|min|plot>

Use the chi2 statistic to choose columns to prune. This will first order
all of the columns by the chi2 statistic by comparing the chi2
statistic for the alignment with and without each column. Then the option
specifies how to remove columns:

    half    Remove half of the sites (starting with the most biased).
    f#      Remove sites until # fraction of sites remains, half can be
            specified as f0.5
    n#      Remove sites until only # number of sequences show significant bias.
    min     Remove sites until a minimum of sequences show significant bias.
    plot    Will print statistics to the screen suitable for plotting. It will
            contain 4 columns: idx, number of biased sequences, chi2 delta for
            this column and the names of the biased sequences.

=cut


has 'chi2_prune' => (
    is            => 'ro',
    isa           => 'Str',
    required      => '',
    default       => '',
    documentation => 'Use the chi2 statistic to prune sites of the alignment, values are Half, n<number>, f<number>, min or Plot.',
);

=head2 --generate_overview

With this option you will get an overview image showing your alignment and the
effect of your pruning settings. The filename will be "infile.png". The colors
are based on BLOSUM62 score spanning from green for really good to red for
really bad, yellow means 0. Removed columns will be overlayed with gray and if
you have specified a subset that will be overlayed with white. Try it to see.

=cut

has 'generate_overview' => (
    is => 'ro',
    isa => 'Bool',
    required => '',
    default => '',
    documentation => 'Generate an overview image of the alignment',
);

=head2 --inverse

Inverse the pruning, that is remove the columns that otherwise would be kept.

=cut

has 'inverse' => (
    is => 'ro',
    isa => 'Bool',
    required => '',
    default => '',
    documentation => 'Inverse the pruning',
);

=head2 --gap_treatment <ignore|mean|additional>

How to treat gaps in the calculations:

    ignore      Just ignore them, default
    mean        Add 1/20 on all the other AAs for each gap in a sequence.
    additional  As an additional state (not implemented)

=cut

has 'gap_treatment' => (
    is => 'ro',
    isa => 'Str',
    required => '',
    default => 'ignore',
    documentation => 'How to treat gaps when calculating conservation and chi2 statistics. Ignore (default), Mean: 1/20 of all other AAs or Additional: additional state',
);

=head2 --aa_significance_level <number>

The significance level used for the chi2 test. The only cases you'd want to
change this is probably when running the chi2 pruning. The default is
30.143527. Can be generated in R with qchisq(0.05, 19, lower.tail=FALSE).

=cut

has aa_significance_level => (
    is      => 'ro',
    isa     => 'Num',
    default => 30.143527,
    documentation => 'The chi2 significance level for amino-acids, default 30.143527',
);

=head2 --dna_significance_level <number>

Same as above but for dna instead. The default is 7.814728. Can be generated in
R with qchisq(0.05, 3, lower.tail=FALSE).

=cut

has dna_significance_level => (
    is      => 'ro',
    isa     => 'Num',
    default => 7.814728,
    documentation => 'The chi2 significance level for nucleotides, default 7.814728',
                              );

=head2 --recode <recoding scheme>

Recoding scheme for recoding if the alignment, currently only dayhoff4 recoding
is supported. Schemes:

    dayhoff4:   A=AGPST  T=DENQ  C=HKR  G=FYWILMV  -=C
    dayhoff6:   1=AGPST  2=DENQ  3=HKR  4=FYW 5=ILMV  6=C
    hp:         1=ACFGILMVW  2=DEHKNPQRSTY

=cut

has recode => (
               is => 'ro',
               isa => 'Str',
               default => '',
               documentation => 'Recoding scheme, dayhoff4',
              );


#### END OF OPTIONS ####

has _matrix => (
    is => 'ro',
    isa => 'HashRef',
    lazy_build => 1,
    traits => ['Hash'],
    handles => {
        _amino_acids => 'keys',
    }
);

sub _build__matrix {
    my $self = shift;

    # The matrix is stored in the DATA-block at the end of this file
    my $matrix = Bio::Matrix::IO->new(
        -fh => \*DATA,
        -format => 'scoring'
    )->next_matrix;

    # Using the matrix directly takes too much time, so we make cache it.
    my %qmatrix;
    my @AAS = grep /[A-Z]/, $matrix->row_names;
    for my $ni ( @AAS ) {
        for my $nj ( @AAS ) {
            $qmatrix{$ni}{$nj} = $matrix->entry($ni,$nj);
        }
    }
    return \%qmatrix;
}


#### MAIN PROGRAM ####

sub run {
    my $self = shift;
    my $alnio = Bio::AlignIO->new(
        -fh => $self->file->openr,
        -format => $self->format,
    );
    $Bio::Root::Root::DEBUG = -1;
    while (my $aln = $alnio->next_aln) {
        $self->prune( $aln );
    }
}

sub prune {
    my $self = shift;
    my ( $aln ) = @_;

    my $stats = $self->calculate_statistics( $aln );
    my $nseq  = $stats->{nseq};

    if ($self->recode) {
        $self->run_recoding( $aln );
        return;
    }

    if ( $self->chi2_test ) {
        $self->run_chi2_test( $stats );
        $self->show_chi2_test( $stats );
        return;
    }

    if ( $self->bowker_symmetry_test ) {
        $self->run_bowker_symmetry_test( $stats );
        return;
    }

    if ( $self->aminogc ) {
        $self->run_aminogc( $stats );
        return;
    }

    if ( $self->chi2_remove_taxa ) {
        $stats = $self->run_chi2_remove_taxa( $aln, $stats );
    }

    my $MAX_COL = $aln->length - 1;
    my @columns_to_remove;

    ## chi2 removal
    if ( $self->chi2_prune ) {
        if ( $self->gap_threshold || $self->conserved_threshold ) {
            die "Can't specify both chi2_prune and any of gap_threshold or conserved_threshold\n";
        }
        # Note that $stats will change after this call, so it's not that usable
        # after this call.
        push @columns_to_remove, $self->run_chi2_prune( $stats );
    }

    ## Gap removal
    if ( $self->gap_threshold ) {
        my $threshold = $self->_convert_threshold( $self->gap_threshold, $nseq );
        for my $pos ( 0 .. $MAX_COL ) {
            next if $stats->{counts}[$pos]{'-'} <= $threshold;
            push @columns_to_remove, $pos;
        }
    }

    ## Unconserved removal
    if ( $self->conserved_threshold ) {
        my $threshold = $self->_convert_threshold( $self->conserved_threshold, $nseq );
        for my $pos ( 0 .. $MAX_COL ) {
            my %counts = %{ $stats->{counts}[$pos] };
            delete $counts{'-'};
            my $best = max values %counts;
            if ( $best < $threshold ) {
                push @columns_to_remove, $pos;
            }
        }
    }

    ## Specific columns removal
    if ( $self->remove_columns ne '' ) {
        my @columns = split /,/, $self->remove_columns;
        for my $c ( @columns ) {
            if ( $c =~ /(\d+)(?:-|\.\.?)(\d+)/ ) {
                push @columns_to_remove, $1 .. $2;
            }
            else {
                push @columns_to_remove, $c;
            }
        }
    }

    if ( $self->inverse ) {
        my %remove = map { ( $_, 1 ) } @columns_to_remove;
        @columns_to_remove =
            grep { ! exists $remove{$_} } 0 .. $aln->length - 1;
    }

    if ( $self->generate_overview ) {
        $self->run_overview_window( $stats, \@columns_to_remove );
    }

#    if ( @columns_to_remove || $self->chi2_remove_taxa ) {
        $aln = $self->run_pruning( $aln, \@columns_to_remove ) if @columns_to_remove;
        $self->write_aln( $aln );
#    }
}

sub write_aln {
    my $self = shift;
    my ($aln) = @_;
    $aln->set_displayname_flat;

    my $IO = Bio::AlignIO->new(
                               -fh => \*STDOUT,
                               -format => $self->format,
                              );

    $IO->write_aln( $aln );
}

sub run_recoding {
    my $self = shift;
    my ($aln) = @_;

    my %schemes = (
                   'dayhoff4' => 'A=AGPST,T=DENQ,C=HKR,G=FYWILMV,-=C',
                   'dayhoff6' => '1=AGPST,2=DENQ,3=HKR,4=FYW,5=ILMV,6=C',
                   'hp'       => '1=ACFGILMVW,2=DEHKNPQRSTY',
                  );
    if (! exists $schemes{ $self->recode }) {
        die "Unrecognized recoding scheme\n";
    }

    my %scheme = (
                  '-' => '-',
                  'X' => 'X',
                 );
    for my $class ( split ',', $schemes{$self->recode}) {
        my ($l,undef,@letters) = split //, $class;
        $scheme{$_} = $l for @letters;
    }

    my @new_seqs;
    for my $seq ( $aln->each_seq ) {
        my $s = $seq->seq;
        $s =~ s{(.)}{$scheme{$1} // die "Invalid alignment character: $1\n"}ge;
        push @new_seqs, Bio::LocatableSeq->new(
                                               -id => $seq->id,
                                               -seq => $s,
                                              );
    }

    $self->write_aln(Bio::SimpleAlign->new(-seqs => \@new_seqs));
}

sub run_pruning {
    my $self = shift;
    my ( $aln, $columns ) = @_;

    return $aln if $columns->[0] == -1;

    my $groups = $self->_group_columns( $columns );
    $aln = $aln->remove_columns( $_ ) for reverse @$groups;

    return $aln;
}

sub run_aminogc {
    my $self = shift;
    my ($stats) = @_;

    my $codon_table = Bio::Tools::CodonTable->new( -id => 11 );
    my @AAS = $self->amino_acids;

    my %chosen;
    for my $AA ( @AAS ) {
        my $codons = join '', map { /^(..)/ } $codon_table->revtranslate( $AA );
        if ( $codons =~ /t|a/i ) {
            $chosen{$AA} = '';
        }
        else {
            $chosen{$AA} = 1;
        }
    }

    my @counts = @{ $stats->{seqcounts} };
    my @names  = @{ $stats->{names} };

    my @aminogc;
    for my $seqi ( 0 .. $#names ) {
        my $aminogc;
        my $tot;

        for my $AA ( keys %{ $counts[$seqi] } ) {
            next if $AA eq '-';
            if ( $chosen{$AA} ) {
                $aminogc += $counts[$seqi]{$AA};
            }
            $tot += $counts[$seqi]{$AA}
        }

        $aminogc /= $tot;
        push @aminogc, $aminogc;

        printf "%-30s %6.4f\n", $names[$seqi], $aminogc;
    }

    say '---';

    my $mean = sum( @aminogc ) / @aminogc;
    my $var  = sum( map { ($_ - $mean)**2 } @aminogc ) / @aminogc;

    printf "%-30s %6.4f\n", 'MEAN', $mean;
    printf "%-30s %6.4f\n", 'STD.DEV.', sqrt($var);
}

sub run_bowker_symmetry_test {
    my $self = shift;
    my ($stats) = @_;

    my $distrib   = $stats->{distribution};
    my $length    = $stats->{length};
    my $nseq      = $stats->{totseq};
    my @AAS       = $self->amino_acids;

    my @names = @{ $stats->{names} };

    my $significance_cutoff = 223.1602;

    my %sign_counts;

    for my $seqi ( 0 .. $nseq - 2 ) {
        for my $seqj ( $seqi+1 .. $nseq - 1 ) {
            my $matrix;
            for my $pos ( 0 .. $length - 1 ) {
                my $aai = $distrib->[$pos][$seqi];
                my $aaj = $distrib->[$pos][$seqj];
                $matrix->{$aai}{$aaj}++
            }
            my $stat = 0;
            for my $ii ( 0 .. $#AAS - 1 ) {
                for my $jj ( $ii+1 .. $#AAS ) {
                    my $ai = $AAS[$ii];
                    my $aj = $AAS[$jj];
                    my $nij = $matrix->{$ai}{$aj} // 0;
                    my $nji = $matrix->{$aj}{$ai} // 0;
                    next unless $nij || $nji;
                    $stat += ( $nij - $nji )**2 / ($nij + $nji);
                }
            }

            printf "%-8.8s %-8.8s %3d %3d %7.2f",
                @names[$seqi,$seqj],
                $seqi, $seqj,
                $stat;
            if ( $stat > $significance_cutoff ) {
                $sign_counts{$names[$seqi]}++;
                $sign_counts{$names[$seqj]}++;
                print ' *';
            }
            print "\n";
        }
    }
}

#### OVERVIEW WINDOW GENERATION ####

sub run_overview_window {
    my $self   = shift;
    my ($stats, $remove) = @_;

    $self->_calculate_statistics_scores( $stats );

    my $distrib   = $stats->{distribution};
    my $scores    = $stats->{scores};
    my $nseq      = $stats->{totseq};
    
    my $AA_WIDTH  = 4;
    my $AA_HEIGHT = 6;

    my $width     = $AA_WIDTH  * $stats->{length};
    my $height    = $AA_HEIGHT * $nseq;

    my $gd = GD::Image->new( $width, $height, 1, );

    ## Cache all colors
    my @COLORS = map {
        $gd->colorAllocate( @$_ )
    } _rainbow_gradient();
    my $white = $gd->colorAllocate( 255,255,255 );

    ## Build the image
    for my $pos ( 0 .. $stats->{length} - 1 ) {
        for my $seqi ( 0 .. $nseq - 1 ) {
            my $aa = $distrib->[$pos][$seqi];
            my $color = $white;
            if ( $aa ne '-' ) {
                my $score = $scores->[$pos][$seqi];
                my $idx = int(100*$score);
                die "Too large $pos,$seqi  $idx,$score" if $idx > 100;
                die "Too small $pos,$seqi  $idx,$score" if $idx < 0;
                $color = $COLORS[ $idx ];
            }

            $gd->filledRectangle(
                $pos      * $AA_WIDTH,
                $seqi     * $AA_HEIGHT,
                ($pos+1)  * $AA_WIDTH,
                ($seqi+1) * $AA_HEIGHT - 1,
                $color,
            );
        }
    }

    ## Add remove overlay
    if ( @$remove ) {
        my $groups = $self->_group_columns( $remove );

        my $color = $gd->colorAllocateAlpha(0,0,0,90);

        for my $group ( @$groups ) {
            my ($start,$end) = @$group;
            $end++;

            $gd->filledRectangle(
                $start * $AA_WIDTH,
                0,
                $end   * $AA_WIDTH,
                $height,
                $color,
            );
        }
    }

    ## Add subset overlay
    if ( $self->subset ) {
        my %chosen = map { ( $_ => 1 ) } @{ $stats->{chosen} };
        my $color = $gd->colorAllocateAlpha(255,255,255,90);        
        for my $c ( 0 .. $nseq - 1 ) {
            next if $chosen{$c};

            $gd->filledRectangle(
                0,
                $c * $AA_HEIGHT,
                $width,
                ($c+1) * $AA_HEIGHT - 1,
                $color,
            );
        }
    }

    my $outfile = $self->file . '.png';
    say STDERR "Saving png to $outfile";
    open my $OUT, '>', $outfile or die;
    print $OUT $gd->png;
    close $OUT;
}


#### CHISQUARE METHODS ####

sub run_chi2_prune {
    my $self = shift;
    my ( $stats ) = @_;

    my $length = $stats->{length};
    my @columns_to_remove;
    given ( $self->chi2_prune ) {
        when ( /^h/i ) {
            my @order = $self->_chi2_prune_order( $stats );
            push @columns_to_remove, @order[ 0 .. $#order * .5];
        }
        when ( /^f(\d+\.?\d*)/i ) {
            my $num   = $self->_convert_threshold( $1, $length );
            my @order = $self->_chi2_prune_order( $stats );
            push @columns_to_remove, @order[ 0 .. ($num-1) ];
        }
        when ( /^n(\d+\.?\d*)/i ) {
            my $num = $self->_convert_threshold( $1, $length );
            push @columns_to_remove,
                $self->_run_chi2_prune_threshold( $stats, $num );
        }
        when ( /^min/i ) {
            push @columns_to_remove,
                $self->_run_chi2_prune_threshold( $stats, 0, 1 );
        }
        when ( /^plot/i ) {
            $self->_plot_chi2_prune( $stats );
            exit;
        }
        default {
            die "chi2 pruning of type $_ not recongnized\n";
        }
    }
    return @columns_to_remove;
}

sub run_chi2_test {
    my $self = shift;
    my ($stats) = @_;

    my @AAS    = $self->amino_acids;
    my @chosen = @{ $stats->{chosen} };
    my @O      = @{ $stats->{seqcounts} };
    my @l      = @{ $stats->{seqlength} };
    my %chosen = map { ($_,1) } @chosen;

    my ( %E );

    for my $seqi ( @chosen ) {
        for my $AA ( @AAS ) {
            $E{$AA} += $O[$seqi]{$AA};
        }
    }

    my $L = sum @l;
    $_ /= $L for values %E;

    my @test;

    my $tot = 0;
    for my $seqi ( 0 .. $#O ) {
        my $chi2 = 0;

        for my $AA ( @AAS ) {
            next if $E{$AA} == 0;
            my $E = $E{$AA} * $l[$seqi];
            $chi2 += ( $O[$seqi]{$AA} - $E )**2 / $E;
        }

        push @test, $chi2;
        $tot += $chi2 if $chosen{$seqi};
    }

    push @test, $tot;
    $stats->{chi2} = \@test;
    return $tot;
}

sub show_chi2_test {
    my $self = shift;
    my ($stats) = @_;

    my $SIGN_LEVEL  = $self->aa_significance_level;

    my @names = @{ $stats->{names} };
    my $chi2  = $stats->{chi2};

    for my $seqi ( 0 .. $#names ) {
        my $sign = $chi2->[$seqi] > $SIGN_LEVEL ? '*' : '';
        
        printf "%4d %-30.30s %7.2f %1s\n",
            $seqi, $names[$seqi], $chi2->[$seqi], $sign;
    }

    printf "     %-30.30s %7.2f\n",
        'TOTAL', $chi2->[-1];

    my $stat = Statistics::Descriptive::Sparse->new;
    $stat->add_data( @$chi2[ 0 .. $#$chi2-1 ] );

    printf "     %-30.30s %7.2f\n",
        'Std.Dev.', $stat->standard_deviation;
}

sub run_chi2_remove_taxa {
    my $self = shift;
    my ($aln, $stats) = @_;

    my $SIGN_LEVEL  = $self->aa_significance_level;

    if ( ! exists $stats->{chi2} ) {
        $self->run_chi2_test( $stats );
    }

    my @names = @{ $stats->{names} };
    my $chi2  = $stats->{chi2};

    my %remove_seqs;

    for my $seqi ( 0 .. $#names ) {
        if ( $chi2->[$seqi] > $SIGN_LEVEL ) {
            $remove_seqs{ $names[$seqi] }++;
        }
    }

    for my $seq ( $aln->each_seq ) {
        if ( exists $remove_seqs{ $seq->id } ) {
            $aln->remove_seq( $seq );
        }
    }
    
    my $new_stats = $self->calculate_statistics( $aln );

    return $new_stats;
}

sub _chi2_test_minus {
    my $self = shift;
    my ($stats, $columns, $inplace) = @_;

    # The inplace flag signifies whether to modify the $stats hash inplace or
    # create a new one.

    my $new_stats = {};

    if ( $inplace ) {
        $new_stats = $stats;
    }
    else {
        $new_stats->{$_} = $stats->{$_} for keys %$stats;
    }

    if ( ref($columns) ne 'ARRAY' ) {
        $columns = [$columns];
    }

    my @distrib  = @{ $stats->{distribution} };
    my @AAS = $self->amino_acids;

    my @O = @{ $stats->{seqcounts} };
    my @l = @{ $stats->{seqlength} };

    my @chosen = @{ $stats->{chosen} };
    my %chosen = map { ($_, 1) } @chosen;

    my (@Oc, @lc, %Ec);

    for my $seqi ( 0 .. $#O ) {
        $lc[$seqi] = $l[$seqi];

        for my $AA ( @AAS ) {
            $Oc[$seqi]{$AA}  = $O[$seqi]{$AA};
            $Ec{$AA}        += $O[$seqi]{$AA} if $chosen{$seqi};
        }

        for my $c ( @$columns ) {
            my $aa = $distrib[$c][$seqi];
            if ( $aa ne '-' ) {
                $Oc[$seqi]{$aa} -= 1;
                $lc[$seqi]      -= 1;
                $Ec{$aa}        -= 1 if $chosen{ $seqi };
            }
        }
    }

    $new_stats->{seqcounts} = \@Oc;
    $new_stats->{seqlength} = \@lc;

    my $L = sum values %Ec;

    if ( $L == 0 ) {
        warn "Nothing left of alignment\n";
        return $new_stats;
    }

    $_ /= $L for values %Ec;

    my @chi2_test;

    my $tot;
    for my $seqi ( @chosen ) {
        my $chi2 = 0;

        for my $AA ( @AAS ) {
            my $E = $Ec{$AA}*$lc[$seqi];
            next if $E == 0;
            $chi2 += ( $Oc[$seqi]{$AA} - $E )**2 / $E;
        }

        $tot += $chi2;
        push @chi2_test, $chi2;
    }
    push @chi2_test, $tot;

    $new_stats->{chi2} = \@chi2_test;
}

sub _chi2_delta_exact {
    my $self = shift;
    my ( $stats, $c, $full_tot ) = @_;

    my @AAS = $self->amino_acids;
    my @distrib  = @{ $stats->{distribution} };

    my @O = @{ $stats->{seqcounts} };
    my @l = @{ $stats->{seqlength} };

    my @chosen = @{ $stats->{chosen} };

    my (@Oc,%Ec,@lc);

    for my $seqi ( @chosen ) {
        $lc[$seqi] = $l[$seqi];

        for my $AA ( @AAS ) {
            $Oc[$seqi]{$AA} = $O[$seqi]{$AA};
            $Ec{$AA} += $O[$seqi]{$AA};
        }

        my $aa = $distrib[$c][$seqi];
        if ( $aa ne '-' ) {
            $Oc[$seqi]{$aa} -= 1;
            $lc[$seqi]      -= 1;
            $Ec{$aa}        -= 1;
        }
    }

    my $L = sum grep defined, @lc;
    $_ /= $L for values %Ec;

    my $tot;
    for my $seqi ( @chosen ) {
        my $chi2 = 0;

        for my $AA ( @AAS ) {
            next if $Ec{$AA} == 0;
            my $E = $Ec{$AA}*$lc[$seqi];
            $chi2 += ( $Oc[$seqi]{$AA} - $E )**2 / $E;
        }

        $tot += $chi2;
    }
    return $tot - $full_tot;
}

sub _chi2_prune_order {
    my $self = shift;
    my ( $stats ) = @_;

    # TODO: Lazy build...
    if ( ! $stats->{chi2} ) {
        $self->run_chi2_test( $stats );
    }

    my $len = $stats->{length};
    my $chi2_sum = $stats->{chi2}[-1];
    my @pos_scores;
    for my $c ( 0 .. $len - 1 ) {
        push @pos_scores, $self->_chi2_delta_exact( $stats, $c, $chi2_sum );
    }

    $stats->{chi2_deltas} = \@pos_scores;

    my @order = sort { $pos_scores[$a] <=> $pos_scores[$b] } 0 .. $len - 1;
    return @order;
}

sub _plot_chi2_prune {
    my $self = shift;
    my ( $stats ) = @_;

    my @order = $self->_chi2_prune_order( $stats );
    @order = sort { $a <=> $b } @order;

    my $new_stats = {};
    $new_stats->{$_} = $stats->{$_} for keys %$stats;
    my $len = $stats->{length};

    my $SIGN_LEVEL = $self->aa_significance_level;
    for my $idx_order ( 0 .. $len - 2 ) {
        $self->_chi2_test_minus( $new_stats, $order[$idx_order], 1 );

        my @delta = @{ $new_stats->{chi2} };
        pop @delta; # Not interested in total

        my @over      = grep { $delta[$_] > $SIGN_LEVEL } 0 .. $#delta;
        my $num       = grep { $_ > $SIGN_LEVEL } @delta;
        my $pos_delta = $stats->{chi2_deltas}[$order[$idx_order]];

        printf "%5d %5d %5.2f (%s)\n", $idx_order, $num, $pos_delta,
            join(',', @{$stats->{names}}[@over]);
    }
}

sub _run_chi2_prune_threshold {
    my $self = shift;
    my ( $stats, $stop_num, $MIN ) = @_;

    $stop_num //= 0;

    my @order = $self->_chi2_prune_order( $stats );

    my $new_stats = {};
    $new_stats->{$_} = $stats->{$_} for keys %$stats;

    my $len = $stats->{length};

    my ($best,$best_idx);

    my $SIGN_LEVEL = $self->aa_significance_level;
    for my $idx_order ( 0 .. $len - 2 ) {
        $self->_chi2_test_minus( $new_stats, $order[$idx_order], 1 );

        my @delta = @{ $new_stats->{chi2} };
        pop @delta; # Not interested in total

        my @over = grep { $delta[$_] > $SIGN_LEVEL } 0 .. $#delta;
        my $num  = grep { $_ > $SIGN_LEVEL } @delta;

        if ( $num <= $stop_num ) {
            return @order[ 0 .. $idx_order ];
        }

        if ( !$best || $num < $best ) {
            ($best,$best_idx) = ($num, $idx_order);
        }
    }

    if ( $MIN ) {
        print STDERR "Min is at $best_idx with $best\n";
        return @order[ 0 .. $best_idx ];
    }

    die "It is not possible to get $stop_num or fewer sequences above significance threshold\n";
}

## Not in current option list
sub _run_reorder_chi2 {
    my $self = shift;
    my ( $stats ) = @_;

    my @order = $self->_chi2_prune_order( $stats );

    my $nseq    = $stats->{nseq};
    my @names   = @{ $stats->{names} };
    my @distrib = @{ $stats->{distribution} };

    my @new_seq_str;
    for my $c ( @order ) {
        for my $seqi ( 0 .. $nseq - 1 ) {
            $new_seq_str[$seqi] .= $distrib[$c][$seqi];
        }
    }

    my @new_seq;
    for my $seqi ( 0 .. $nseq - 1 ) {
        push @new_seq, Bio::LocatableSeq->new(
            -id => $names[$seqi],
            -seq => $new_seq_str[$seqi],
        );
    }

    return Bio::SimpleAlign->new(
        -seqs => \@new_seq,
    );
}


# TODO: Should be replaced with a different sub to return the alphabet of the
#       alignment instead, that way we can support dna quite easily.

sub amino_acids {
    my $self = shift;
    if ( $self->gap_treatment =~ /^a/i ) {
        return ( $self->_amino_acids, '-');
    }
    return $self->_amino_acids;
}

#### STATISTICS ####

sub calculate_statistics {
    my $self = shift;
    my $aln  = shift;

    # TODO: This should be refactored into it's own class, but I don't have the
    #       tuits...
    #
    # In the end we have the following datastructure. 
    #$stats = {
    #### Indexed by column
    #    scores       Array of Arrays [pos][seqi]  = Normalized Blosum62 score
    #    distribution Array of Arrays [pos][seqi]  = AA at pos in seqi
    #    counts       Array of Hashes [pos]{char}  = n(AA) in column pos
    #
    #### Indexed by sequence
    #    seqcounts    Array of Hashes [seqi]{char} = n(AA) in seqi
    #    frequency    Array of Hashes [seqi]{char} = f(AA) in column sequence seqi, between 0 and 1
    #    names        Array           [seqi]       = name of seqi
    #    seqlength    Array           [seqi]       = l(seqi)
    #
    #### Other stuff
    #    chosen       Array           [idx]        = chosen seqi
    #    nseq         Integer                      = number of chosen seqs
    #    totseq       Integer                      = total number of seqs
    #    length       Integer                      = length of alignment
    #    overall_freq Hash            {char}       = Overall frequency
    #};


    # Each of these methods builds upon the previous one.

    my $stats = {};
    $self->_calculate_statistics_transpose( $stats, $aln );
    $self->_calculate_statistics_counts( $stats );
    # This one takes a lot of time, we defer it until it is needed
    #$self->_calculate_statistics_scores( $stats );

    return $stats;
}

sub _calculate_statistics_transpose {
    my $self = shift;
    my ($stats, $aln) = @_;

    my $subset = $self->subset;
    my $totseq = $aln->num_sequences;
    my @AAS = $self->amino_acids;

    my (@distrib, @names, @chosen, @seqcounts);
    my (@frequency, %overall_frequency, @length);

    my $seqi=0;
    for my $seq ( $aln->each_seq ) {
        my $chosen = '';
        if ( !$subset || $seq->id =~ /$subset/ ) {
            $chosen++;
            push @chosen, $seqi;
        }

        push @names, $seq->id;

        my %seqcounts = map { ($_=>0) } ( @AAS, '-' );
        my @s = split //, $seq->seq;
        for my $pos ( 0 .. $#s ) {
            $seqcounts{$s[$pos]}++;
            push @{ $distrib[$pos] }, $s[$pos];
        }

        given ( $self->gap_treatment ) {
            when (/^i/i) { delete $seqcounts{'-'} }
            when (/^m/i) {
                $seqcounts{$_} += $seqcounts{'-'}/20 for @AAS;
                delete $seqcounts{'-'};
            }
            when (/^a/i) {
                die "Treatment of gaps as an additional character is not implemented yet\n";
            }
        }

        push @seqcounts, \%seqcounts;

        if ( $chosen ) {
            $overall_frequency{$_} += $seqcounts{$_} for keys %seqcounts;
        }
        
        my $sum = sum values %seqcounts;
        push @length, $sum;

        my %seq_frequency = %seqcounts;
        $_ /= $sum for values %seq_frequency;
        push @frequency, \%seq_frequency;

        $seqi++;
    }

    my $sum = sum values %overall_frequency;
    $_ /= $sum for values %overall_frequency;

    $stats->{distribution} = \@distrib;
    $stats->{chosen}       = \@chosen;
    $stats->{names}        = \@names;
    $stats->{totseq}       = $totseq;
    $stats->{nseq}         = scalar(@chosen);
    $stats->{length}       = scalar(@distrib);
    $stats->{seqcounts}    = \@seqcounts;
    $stats->{seqlength}    = \@length;
    $stats->{frequency}    = \@frequency;
    $stats->{overall_freq} = \%overall_frequency;
}

sub _calculate_statistics_counts {
    my $self = shift;
    my ($stats) = @_;

    my @AAS = ($self->amino_acids, '-');
    my @chosen = @{ $stats->{chosen} };

    my (@counts);

    for my $p ( @{ $stats->{distribution} } ) {
        my %counts = map { ($_ => 0) } @AAS;

        $counts{$_}++ for @$p[@chosen];
        push @counts, \%counts;
    }

    $stats->{counts} = \@counts;
}

sub _calculate_statistics_scores {
    my $self = shift;
    my ( $stats ) = @_;

    my $qmatrix = $self->_matrix;
    my @AAS = $self->amino_acids;

    # Precalculate the scores as weights
    my @weights;
    for my $pos ( 0 .. $stats->{length} - 1 ) {
        my $p = $stats->{distribution}[$pos];

        my %weights;

        for my $aai ( @AAS ) {
            my $score = 0;

            for my $aaj ( @AAS ) {
                if ( ! defined $stats->{counts}[$pos]{$aaj} ) {
                    die "Can't find $aaj at $pos of counts";
                }
                $score += $qmatrix->{$aai}{$aaj} * $stats->{counts}[$pos]{$aaj};
            }

            # Remove self comparison
            $score -= $qmatrix->{$aai}{$aai};

            $weights{$aai} = $score;
        }
        my ($smin, $smax) = minmax values %weights;

        # Normalize the scores so they are mapped like this:
        #     (-Inf,0,Inf) -> (0 , 0.5 , 1)
        for ( values %weights ) {
            when ( $_    < 0 ) { $_ = (1 - $_/$smin) / 2 }
            when ( $smax > 0 ) { $_ = 0.5 + ($_/$smax) / 2 }
        }

        push @weights, \%weights;
    }

    # Map the weights to the scores
    my @scores;
    for my $pos ( 0 .. $#weights ) {
        for my $seqi ( 0 .. $stats->{totseq}-1 ) {
            my $c = $stats->{distribution}[$pos][$seqi];
            next if $c eq '-';

            $scores[$pos][$seqi] = $weights[$pos]{$c};
        }
    }

    $stats->{scores} = \@scores;
}

#### UTILITY FUNCTIONS ####

sub _group_columns {
    my $self = shift;
    my $columns = shift;

    my @columns = sort { $a <=> $b } uniq @$columns;
    my @groups;
    my $current_group = [];
    for ( @columns ) {
        if ( !@$current_group ) {
            $current_group = [$_,$_];
        }
        elsif ( $current_group->[1] + 1 == $_ ) {
            $current_group->[1] = $_;
        }
        else {
            push @groups, $current_group;
            $current_group = [$_,$_]
        }
    }
    push @groups, $current_group;
    return \@groups;
}

sub _convert_threshold {
    my $self = shift;
    my ( $value, $nseq ) = @_;

    my $ret;
    for ( $value ) {
        when (/^(\d+\.?\d*)\%$/) { $ret = $nseq * $1 / 100 }
        when (/^(\d+)$/)         { $ret = $value }
        when (/^(\d*\.\d+)$/)    { $ret = $nseq * $1 }
        default { die "Can't interpret threshold: $value\n" }
    }

    if ( $ret > $nseq ) {
        die "Threshold can't be larger than nseq $value ($ret) > $nseq\n";
    }

    return $ret;
}

sub _gradient {
    map { [ map hex, /(..)(..)(..)/ ] } @_;
}

sub _jalview_gradient {
    # Generated in R with hsv(2/3, seq(1,101)/101, 1)
    _gradient qw(
        FCFCFF FAFAFF F7F7FF F5F5FF F2F2FF F0F0FF EDEDFF EBEBFF E8E8FF E6E6FF
        E3E3FF E1E1FF DEDEFF DCDCFF D9D9FF D7D7FF D4D4FF D2D2FF CFCFFF CDCDFF
        CACAFF C7C7FF C5C5FF C2C2FF C0C0FF BDBDFF BBBBFF B8B8FF B6B6FF B3B3FF
        B1B1FF AEAEFF ACACFF A9A9FF A7A7FF A4A4FF A2A2FF 9F9FFF 9D9DFF 9A9AFF
        9797FF 9595FF 9292FF 9090FF 8D8DFF 8B8BFF 8888FF 8686FF 8383FF 8181FF
        7E7EFF 7C7CFF 7979FF 7777FF 7474FF 7272FF 6F6FFF 6D6DFF 6A6AFF 6868FF
        6565FF 6262FF 6060FF 5D5DFF 5B5BFF 5858FF 5656FF 5353FF 5151FF 4E4EFF
        4C4CFF 4949FF 4747FF 4444FF 4242FF 3F3FFF 3D3DFF 3A3AFF 3838FF 3535FF
        3232FF 3030FF 2D2DFF 2B2BFF 2828FF 2626FF 2323FF 2121FF 1E1EFF 1C1CFF
        1919FF 1717FF 1414FF 1212FF 0F0FFF 0D0DFF 0A0AFF 0808FF 0505FF 0303FF
        0000FF
    )
}

sub _rainbow_gradient {
    # Generated in R with rainbow(101, s=.6, v=.9, start=0, end=1/3)
    _gradient qw(
        E65C5C E65F5C E6615C E6645C E6675C E66A5C E66C5C E66F5C E6725C E6755C
        E6775C E67A5C E67D5C E6805C E6825C E6855C E6885C E68B5C E68D5C E6905C
        E6935C E6965C E6985C E69B5C E69E5C E6A15C E6A35C E6A65C E6A95C E6AC5C
        E6AE5C E6B15C E6B45C E6B75C E6B95C E6BC5C E6BF5C E6C25C E6C45C E6C75C
        E6CA5C E6CD5C E6CF5C E6D25C E6D55C E6D85C E6DA5C E6DD5C E6E05C E6E35C
        E6E65C E3E65C E0E65C DDE65C DAE65C D8E65C D5E65C D2E65C CFE65C CDE65C
        CAE65C C7E65C C4E65C C2E65C BFE65C BCE65C B9E65C B7E65C B4E65C B1E65C
        AEE65C ACE65C A9E65C A6E65C A3E65C A1E65C 9EE65C 9BE65C 98E65C 96E65C
        93E65C 90E65C 8DE65C 8BE65C 88E65C 85E65C 82E65C 80E65C 7DE65C 7AE65C
        77E65C 75E65C 72E65C 6FE65C 6CE65C 6AE65C 67E65C 64E65C 61E65C 5FE65C
        5CE65C
    )
}


__PACKAGE__->new_with_options->run unless caller;


### Don't change anything below this line unless you know how it's used
__DATA__
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1 -1 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -1 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -1 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0 -1 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1 -1 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -1 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
