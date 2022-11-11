#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Term::ANSIColor 4.00; #Needs >= 4.0 to use 256 colors
use Carp qw(confess); #use confess instead of die to print the traceback of error, useful sometimes
use Pod::Usage;
use Getopt::Long;


my $help;
my $input;
my $phredscore = 33;
my $bargraph = 1; #Default output
my $score;
my $blackwhite;
my $html;
my $maxseqs;

GetOptions(
    'h|help!' => \$help,
    'i|input=s' => \$input,
    'p|phredscore=i' => \$phredscore,
    'b|bargraph' => \$bargraph,
    's|score' => \$score,
    'w|blackwhite' => \$blackwhite,
    'o|out=s' => \$html,
    'm|maxseqs=i' => \$maxseqs,
);

pod2usage( -verbose => 2) if ( $help );

unless ( $input ) {
    pod2usage(
    	-msg => "** ERROR: Invalid input file $! **",
    	-verbose => 1
    );
}
unless ($phredscore == 33 || $phredscore == 64) {
    pod2usage(
    	-msg => "** ERROR: Phredscore can only be '33' or '64' **",
    	-verbose => 1
	);
}

my $in  = Bio::SeqIO->new(-format => 'fastq',
			  -file   => $input);

#Check if a HTML output is asked for. If not, print to STDOUT
if ($html) {
    #Load extra modules needed
    use Color::ANSI::Util qw(ansi256_to_rgb);
    open(OUT, '>:utf8', $html) or die;
    print OUT &htmlstart($input); #Prints starting blocks of HTML file
} else {
    open(OUT, '>&:utf8', \*STDOUT) or die; #utf8 so colors and characters work
}

#Color phredscore
if ($score) {
    my $seqcounter = 0;
    while (my $seq = $in->next_seq) {
	#Print the Header
	print OUT "@" . $seq->id;
	print OUT " " . $seq->desc if ($seq->desc);
	print OUT $html ? "<br>\n" : "\n";
	#Print the sequence
	print OUT $seq->seq;
	print OUT $html ? "<br>\n" : "\n";
	#Print a '+' as quality header
	print OUT $html ? "+<br>\n" : "+\n";
	my @arr = split(/ /,$seq->qual_text);
	chomp(@arr);
	my $xcount = 1;
	foreach my $x (@arr) {
	    print STDERR $xcount, "\t", $x, "\n";
	    $xcount++;
	}
	for (@arr) {
	    if ($phredscore == 33) {
		if ($html) {
		    my ($red, $green, $blue) = (split(//, &phred2col($_)))[3,4,5];
		    #Convert to proper RGB value
		    my $rgb = 16 + 36 * $red + 6 * $green + $blue;
		    #Get hex value
		    my $hex = ansi256_to_rgb($rgb);
		    print OUT "<font color=\"#$hex\">";
		    print OUT chr($_+33);
		    print OUT "</font>";
		}
		else {
		    print OUT color(&phred2col($_)) unless ($blackwhite || $html);
		    print OUT chr($_+33);
		}
	    }
	    elsif ($phredscore == 64) {
		if ($html) {
		    my ($red, $green, $blue) = (split(//, &phred2col($_-33)))[3,4,5];
		    #Convert to proper RGB value
		    my $rgb = 16 + 36 * $red + 6 * $green + $blue;
		    #Get hex value
		    my $hex = ansi256_to_rgb($rgb);
		    print OUT "<font color=\"#$hex\">";
		    print OUT chr($_+33);
		    print OUT "</font>";
		}
		else {
		    print OUT color(&phred2col($_-33)) unless ($blackwhite || $html);
		    print OUT chr($_+33);
		}
		print OUT color(&phred2col($_-33)) unless ($blackwhite || $html);
		print OUT chr($_+33);
	    }
	}
	print OUT color ('reset') unless ($blackwhite || $html);
	print OUT $html ? "<br>\n" : "\n";
	#print OUT "\n";
	$seqcounter++;
	if ($maxseqs && $seqcounter > $maxseqs) {
	    last;
	}
    }
}

#Bar-diagram
elsif ($bargraph) {
    my $seqcounter = 0;
    while (my $seq = $in->next_seq) {
	#Print the Header
	print OUT "@" . $seq->id;
	print OUT " " . $seq->desc if ($seq->desc);
	print OUT $html ? "<br>\n" : "\n";
	#Print the sequence
	print OUT $seq->seq;
	print OUT $html ? "<br>\n" : "\n";
	#Print a '+' as quality header
	print OUT $html ? "+<br>\n" : "+\n";
	#print OUT "+<br>\n" if $html;
	my @arr = split(/ /,$seq->qual_text);
	for (@arr) {
	    if ($phredscore == 33) {
		if ($html) {
		    my ($red, $green, $blue) = (split(//, &phred2barsCol($_)))[3,4,5];
		    #Convert to proper RGB value
		    my $rgb = 16 + 36 * $red + 6 * $green + $blue;
		    #Get hex value
		    my $hex = ansi256_to_rgb($rgb);
		    print OUT "<font color=\"#$hex\">";
		    print OUT &phred2bars($_);
		    print OUT "</font>";
		}
		else {
		print OUT color(&phred2barsCol($_)) unless ($blackwhite || $html);
		print OUT &phred2bars($_);
		}
	    }
	    elsif ($phredscore == 64) {
		if ($html) {
		    my ($red, $green, $blue) = (split(//, &phred2barsCol($_-33)))[3,4,5];
		    #Convert to proper RGB value
		    my $rgb = 16 + 36 * $red + 6 * $green + $blue;
		    #Get hex value
		    my $hex = ansi256_to_rgb($rgb);
		    print OUT "<font color=\"#$hex\">";
		    print OUT &phred2bars($_-33);
		    print OUT "</font>";
		}
		else {
     		print OUT color(&phred2barsCol($_-33)) unless ($blackwhite || $html);
		print OUT &phred2bars($_-33);
		}
	    }	
	}
	print OUT color ('reset') unless ($blackwhite || $html);
	print OUT $html ? "<br>\n" : "\n";

	$seqcounter++;
	if ($maxseqs && $seqcounter > $maxseqs) {
	    last;
	}
    }
}

close OUT;

#Really ugly subs start here

sub phred2bars {
    my ($nScore) = @_;
    if ($nScore < 3) {
	return '_';
    }
    elsif ($nScore >= 3 && $nScore < 7) {
	return "\x{2581}";
    }
    elsif ($nScore >= 7 && $nScore < 11) {
	return "\x{2582}";
    }
    elsif ($nScore >= 11 && $nScore < 15) {
	return "\x{2583}";
    }
    elsif ($nScore >= 15 && $nScore < 19) {
	return "\x{2584}";
    }
    elsif ($nScore >= 19 && $nScore < 23) {
	return "\x{2585}";
    }
    elsif ($nScore >= 23 && $nScore < 27) {
	return "\x{2586}";
    }
   elsif ($nScore >= 27 && $nScore < 31) {
	return "\x{2587}";
    }
   elsif ($nScore >= 31) {
	return "\x{2588}";
    }
    else {
	die "Non valid phredscore!\n";
    }
}

sub phred2barsCol {
    my ($nScore) = @_;
    if ($nScore < 3) {
	return "RGB100";
    }
    elsif ($nScore >= 3 && $nScore < 7) {
	return "RGB510";
    }
    elsif ($nScore >= 7 && $nScore < 11) {
	return "RGB520";
    }
    elsif ($nScore >= 11 && $nScore < 15) {
	return "RGB530";
    }
    elsif ($nScore >= 15 && $nScore < 19) {
	return "RGB540";
    }
    elsif ($nScore >= 19 && $nScore < 23) {
	return "RGB440";
    }
    elsif ($nScore >= 23 && $nScore < 27) {
	return "RGB340";
    }
    elsif ($nScore >= 27 && $nScore < 31) {
	return "RGB240";
    }
    elsif ($nScore >= 31) {
	return "RGB040";
    }
    else {
	die "Non valid phredscore!\n";
    }
}

sub phred2col {
    my ($nScore) = @_;
    if ($nScore < 3) {
	return "RGB400";
    }
    elsif ($nScore >= 3 && $nScore < 6) {
	return "RGB510";
    }
    elsif ($nScore >= 6 && $nScore < 9) {
	return "RGB520";
    }
    elsif ($nScore >= 9 && $nScore < 12) {
	return "RGB530";
    }
    elsif ($nScore >= 12 && $nScore < 15) {
	return "RGB540";
    }
    elsif ($nScore >= 15 && $nScore < 18) {
	return "RGB440";
    }
    elsif ($nScore >= 18 && $nScore < 21) {
	return "RGB340";
    }
   elsif ($nScore >= 21 && $nScore < 24) {
	return "RGB240";
    }
   elsif ($nScore >= 24 && $nScore < 27) {
	return "RGB140";
    }
   elsif ($nScore >= 27 && $nScore < 30) {
	return "RGB040";
    }
   elsif ($nScore >= 30) {
	return "RGB050";
    }
    else {
	die "Non valid phredscore!\n";
    }
}

sub htmlstart {
    my ($inputname) = @_;
    return "<HTML>\n<HEAD>\n<TITLE>$inputname</TITLE>\n</HEAD>\n<style type=\"text/css\">\nbody{\nwhite-space: nowrap;\n}\n</style>\n<BODY bgcolor=\"#000000\"><meta http-equiv=\"Content-Type\" content=\"text/html;charset=UTF-8\">\n<FONT COLOR=\"#FFFFFF\" FACE=\"Courier New\">\n";
}

sub htmlend {
    return "</FONT></BODY>\n</HTML>\n";
}

=head1 NAME

    colorFastq.pl - Outputs fastq with colorized phred-scores or bargraphs.

=head1 USAGE

    colorFastq.pl -i input-fastq -b (output bargraph [default] -s (color score only) -w (black/white)

=head1 DESCRIPTION

    Outputs colorized fastq with bargraphs or 'raw' ascii values. 
    Useful for getting a better understanding of the quality of individual reads.
    Will print results to STDOUT.

=head1 OPTIONS

=head2 REQUIRED

=over 3

=item -i, --input
    input [file], fastq format

=back

=head2 OPTIONAL

=over 3

=item -p, --phredscore 
    Input phredscore type
    Can be '33' or '64' (default: phred33)
    
=item -b, --bargraph 
    Output bargraphs (default)

=item -s, --score 
    Output colored phredscores

=item -w,
    Disable the colors. 
    Useful if used on a termianl that does not support colors
    If used together with -s, will produce a normal looking fastq file
    Why would you even use this script then?

=item -h, --help
    Shows the perldoc

=back

=cut

=head1 DEPENDENCIES

    Term::ANSIColor, version > 4.0
   
=head1 AUTHOR

    Anders Lind, anders.lind@icm.uu.se

=cut
