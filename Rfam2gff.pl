#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

=head1 NAME

Rfam2gff.pl - Convert Rfam search output to GFF3 format

=head1 SYNOPSIS

    perl Rfam2gff.pl [options] < rfam_output.txt > output.gff3
    
Options:
    --evalue    E-value threshold for filtering (default: 0.001)
    --help      Display this help message
    --man       Display detailed documentation

=head1 DESCRIPTION

This script converts Rfam search output to GFF3 format. It processes the input
line by line, filtering by E-value, and generates GFF3-compliant output with
gene and exon features. Each feature is given a unique ID based on sequence
name and hit number.

=cut

# Default parameters
my $evalue_threshold = 0.001;
my $help = 0;
my $man = 0;

# Parse command line options
GetOptions(
    'evalue=f' => \$evalue_threshold,
    'help'     => \$help,
    'man'      => \$man
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;

# Track number of occurrences for each model
my %model_counts;

# Print GFF3 header
print "##gff-version 3\n";

while (my $line = <STDIN>) {
    chomp $line;
    
    # Skip comment lines
    next if $line =~ /^#/;
    
    # Parse line into fields
    my @fields = split(/\s+/, $line);
    
    # Validate input format
    unless (@fields >= 16) {
        warn "Warning: Skipping malformed line: $line\n";
        next;
    }
    
    # Extract relevant fields
    my ($seq_id, $model_id, $start, $end, $score, $e_value) = 
        @fields[0, 2, 7, 8, 14, 15];
    
    # Filter by E-value threshold
    next if $e_value > $evalue_threshold;
    
    # Track model occurrences
    $model_counts{$model_id} = 0 unless exists $model_counts{$model_id};
    $model_counts{$model_id}++;
    
    # Determine strand and ensure start < end
    my ($final_start, $final_end, $strand);
    if ($start <= $end) {
        $final_start = $start;
        $final_end = $end;
        $strand = '+';
    } else {
        $final_start = $end;
        $final_end = $start;
        $strand = '-';
    }
    
    # Generate feature IDs
    my $feature_base = "$seq_id.$model_id." . $model_counts{$model_id};
    my $gene_id = "TU.$feature_base";
    my $exon_id = "exon.$feature_base";
    
    # Print GFF3 features
    print join("\t",
        $seq_id, 'Rfam', 'gene',
        $final_start, $final_end, $score,
        $strand, '.',
        "ID=$gene_id;Name=$gene_id"
    ) . "\n";
    
    print join("\t",
        $seq_id, 'Rfam', 'exon',
        $final_start, $final_end, $score,
        $strand, '.',
        "ID=$exon_id;Parent=$gene_id"
    ) . "\n";
}
