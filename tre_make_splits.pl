#!/usr/bin/env perl
package make_splits;
use strict;
use warnings;
use feature ':5.10';
use Moose;
use MooseX::Types::Path::Class;
with 'MooseX::Getopt';

use Phylogeny::Tree::Split;
use Phylogeny::Tree::SplitSet;
use Phylogeny::Tree::Parser;
use Phylogeny::Tree::Node;

has 'file' => (
    is => 'ro',
    isa => 'Path::Class::File',
    required => 1,
    coerce => 1,
);

has 'outfile' => (
    is => 'ro',
    isa => 'Path::Class::File',
    coerce => 1,
    required => 1,
);

sub run {
    my $self = shift;

    my $OUT = $self->outfile->openw;
    my $split_set = $self->get_split_set();
    my $str = $split_set->stringify;
    #$str =~ s/\|[0-9A-F]+//ig;
    say $OUT $str;
}

sub get_split_set {
    my $self = shift;

    my $IO_TREE = $self->file->openr;

    my $split_set = Phylogeny::Tree::SplitSet->new;
    while ( my $line = <$IO_TREE> ) {
        $line =~ s/;\s*$//g;
        my $tree = Phylogeny::Tree::Parser->parse($line);
        my @splits = $tree->get_splits;
        $split_set->add_split($_) for @splits;
    }

    return $split_set;
}

make_splits->new_with_options->run unless caller;

1;
