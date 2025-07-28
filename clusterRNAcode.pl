#!/usr/bin/perl
# ===========================================================================
# clusterRNAcode.pl
# ===========================================================================
#
# Description:
#   Clusters and merges overlapping RNAcode predictions from sorted input.
#   Takes RNAcode output sorted by chromosome, start, and end positions,
#   and merges overlapping regions while tracking best scores and p-values.
#
# Usage:
#   cat rnaCode.xx* | sort -k 1,1 -k 4,4g -k 5,5g | grep Supercontig | perl clusterRNAcode.pl
#
# Input Format (tab-delimited):
#   chr source feature start end score|pvalue strand frame attributes
#
# Output Format (GTF):
#   chr RNAcode CDS start end score|pvalue strand . gene_id "GeneX"
#
# Dependencies:
#   Perl 5.10 or higher
#   Data::Dumper (for debugging)
#
# Author: htafer
# Last Updated: 2025-07-28
# ===========================================================================

use strict;
use warnings;
use feature 'say';
use Data::Dumper;

# Constants
use constant {
    MIN_SCORE => 1,
    MIN_PVALUE => 1,
    DEFAULT_STRAND => "/",
    EMPTY_CHR => "",
    INVALID_POS => -1,
};

# Initialize locus data structure with default values
my $locus = {
    prChr => EMPTY_CHR,
    prStart => INVALID_POS,
    prStop => INVALID_POS,
    prStrand => DEFAULT_STRAND,
    prScore => MIN_SCORE,
    prPValue => MIN_PVALUE,
    count => 0,
};
# Subroutine to print a locus in GTF format
sub print_locus {
    my ($locus) = @_;
    return if $locus->{prChr} eq EMPTY_CHR;  # Skip empty initial locus
    
    say join("\t",
        $locus->{prChr},
        "RNAcode",
        "CDS",
        $locus->{prStart},
        $locus->{prStop},
        $locus->{prScore} . "|" . $locus->{prPValue},
        $locus->{prStrand},
        ".",
        "gene_id \"Gene" . $locus->{count} . "\""
    );
}

# Process input line by line
while (my $line = <STDIN>) {
    chomp $line;
    
    # Parse input line
    my @data = split(/\t/, $line);
    die "Invalid input format: expected 9 fields, got " . scalar(@data) 
        unless @data >= 7;
    
    # Extract score and p-value
    my ($pValue, $score) = split(/\|/, $data[5]);
    die "Invalid score format: $data[5]" unless defined $pValue && defined $score;
    
    # Check if current feature should start a new locus
    if ($data[0] ne $locus->{prChr} ||
        $data[3] > $locus->{prStop} ||
        $data[6] ne $locus->{prStrand}) {
        
        # Print previous locus
        print_locus($locus);
        $locus->{count}++;
        
        # Initialize new locus
        $locus->{prChr} = $data[0];
        $locus->{prStart} = int($data[3]);
        $locus->{prStop} = int($data[4]);
        $locus->{prStrand} = $data[6];
        $locus->{prScore} = $score;
        $locus->{prPValue} = $pValue;
    }
    else {
        # Merge with existing locus
        $locus->{prStart} = int($data[3]) if int($data[3]) < $locus->{prStart};
        $locus->{prStop} = int($data[4]) if int($data[4]) > $locus->{prStop};
        $locus->{prPValue} = $pValue if $pValue < $locus->{prPValue};
        $locus->{prScore} = $score if $score < $locus->{prScore};
    }
}

# Print final locus
print_locus($locus);


