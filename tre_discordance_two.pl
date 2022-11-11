#!/usr/bin/env perl
package tre_discordance_two;
use strict;
use warnings;
use feature ':5.10';
use Moose;
use MooseX::Types::Path::Class;
with 'MooseX::Getopt';

use FindBin;
use lib "$FindBin::Bin";
use Phylogeny::Tree::SplitSet;

has 'file1' => (
    is => 'ro',
    isa => 'Path::Class::File',
    coerce => 1,
    required => 1,
    documentation => 'First file of splits',
);

has 'file2' => (
    is => 'ro',
    isa => 'Path::Class::File',
    coerce => 1,
    required => 1,
    documentation => 'Second file of splits',
);

has outfile => (
    is => 'ro',
    isa => 'Path::Class::File',
    coerce => 1,
    required => '',
    documentation => 'Outfile to put result in',
);

has cutoff => (
    is            => 'ro',
    isa           => 'Int',
    required      => '',
    default       => 75,
    documentation => 'Cutoff to use for discordance filter (default 75)',
);


sub run {
    my $self = shift;

    my $split_set1 = $self->get_split_set( $self->file1 );
    my $split_set2 = $self->get_split_set( $self->file2 );
    my $conflict_score = $split_set1->conflict_score( $split_set2, $self->cutoff, '' );

    my $OUT = $self->outfile->openw;
    say $OUT $conflict_score;
}

sub get_split_set {
    my $self = shift;
    my $splits_file = shift;

    open my $IO_TREE, '<', $splits_file or die "Can't open file $splits_file: $!";
    my $str = join '', <$IO_TREE>;

    my $split;
    eval {
        $split = Phylogeny::Tree::SplitSet->from_string( $str );
    }; if ( $@ ) {
        die "Can't parse split from $splits_file: $@";
    }
    return $split;
}

tre_discordance_two->new_with_options->run unless caller;

__PACKAGE__->meta->make_immutable;

1;
