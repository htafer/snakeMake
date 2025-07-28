#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

=head1 NAME

tRNAscan2gff.pl - Convert tRNAscan-SE output to GFF3 format

=head1 SYNOPSIS

    perl tRNAscan2gff.pl [options] < tRNAscan-SE.out > output.gff3
    
Options:
    --skip-header    Number of header lines to skip (default: 3)
    --help          Display this help message
    --man           Display detailed documentation

=head1 DESCRIPTION

This script converts tRNAscan-SE output to GFF3 format. It processes the input
line by line and generates GFF3-compliant output with gene and exon features.
It handles both single-exon and split tRNAs (intron-containing) correctly.

Input format expected from tRNAscan-SE:
Column 1: Sequence name
Column 2: tRNA number
Column 3: Start position
Column 4: End position
Column 5: tRNA type
Column 6,7: Intron coordinates (if present)
Column 8: Score

=cut

# Default parameters
my $header_lines = 3;
my $help = 0;
my $man = 0;

# Parse command line options
GetOptions(
    'skip-header=i' => \$header_lines,
    'help'         => \$help,
    'man'          => \$man
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;

# Print GFF3 header
print "##gff-version 3\n";

# Skip header lines
for (1..$header_lines) {
    my $header = <STDIN>;
    unless (defined $header) {
        die "Error: Input file appears to be empty or too short\n";
    }
}

sub print_gff_feature {
    my ($seq_id, $source, $type, $start, $end, $score, $strand, $phase, $attrs) = @_;
    print join("\t",
        $seq_id, $source, $type, $start, $end, $score, $strand, $phase, $attrs
    ) . "\n";
}

while (my $line = <STDIN>) {
    chomp $line;
    
    # Parse line into fields
    my @fields = split(/\s+/, $line);
    
    # Validate input format
    unless (@fields >= 9) {
        warn "Warning: Skipping malformed line: $line\n";
        next;
    }
    
    # Extract relevant fields
    my ($seq_id, $trna_num, $start, $end, $trna_type, $intron_start, $intron_end, undef, $score) = @fields;
    
    # Determine strand and ensure start < end
    my ($final_start, $final_end, $strand);
    if ($start < $end) {
        $final_start = $start;
        $final_end = $end;
        $strand = '+';
    } else {
        $final_start = $end;
        $final_end = $start;
        $strand = '-';
    }
    
    # Generate feature ID
    my $feature_base = "$seq_id.tRNA-$trna_type.$trna_num";
    my $gene_id = "TU.$feature_base";
    
    # Print gene feature
    print_gff_feature(
        $seq_id, 'tRNAscan-SE', 'gene',
        $final_start, $final_end, $score, $strand, '.',
        "ID=$gene_id;Name=$gene_id"
    );
    
    # Handle exons
    if ($intron_start + $intron_end > 0) {
        # Split tRNA with intron
        my ($exon1_start, $exon1_end, $exon2_start, $exon2_end);
        if ($strand eq '+') {
            ($exon1_start, $exon1_end) = ($final_start, $intron_start);
            ($exon2_start, $exon2_end) = ($intron_end, $final_end);
        } else {
            ($exon1_start, $exon1_end) = ($intron_end, $final_end);
            ($exon2_start, $exon2_end) = ($final_start, $intron_start);
        }
        
        # Print split exons
        print_gff_feature(
            $seq_id, 'tRNAscan-SE', 'exon',
            $exon1_start, $exon1_end, $score, $strand, '.',
            "ID=exon.$feature_base.exon1;Parent=$gene_id"
        );
        print_gff_feature(
            $seq_id, 'tRNAscan-SE', 'exon',
            $exon2_start, $exon2_end, $score, $strand, '.',
            "ID=exon.$feature_base.exon2;Parent=$gene_id"
        );
    } else {
        # Single exon tRNA
        print_gff_feature(
            $seq_id, 'tRNAscan-SE', 'exon',
            $final_start, $final_end, $score, $strand, '.',
            "ID=exon.$feature_base.exon1;Parent=$gene_id"
        );
    }
}
