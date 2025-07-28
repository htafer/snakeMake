#!/usr/bin/perl
# ====================================
# revigo.pl - REVIGO Web Interface CLI
# ====================================
#
# Purpose:
#   Command-line interface to the REVIGO (REduce + VIsualize Gene Ontology) web service
#   http://revigo.irb.hr/
#
# Description:
#   This script automates the submission of GO terms and p-values to REVIGO,
#   which reduces redundancy in GO term lists and provides visualizations.
#
# Usage:
#   revigo.pl -f <input_file> [options]
#
# Options:
#   -f <file>   Input file containing GO terms and p-values
#   -c <float>  Cutoff for clustering (default: 0.7)
#   -i <yes|no> Second column is p-value (default: yes)
#   -w <string> Value ordering (higher/lower, default: higher)
#   -g <int>    GO universe size (default: 4932 for yeast)
#   -m <string> Similarity measure (default: SIMREL)
#   -t <string> GO type (BP|CC|MF|all, default: BP)
#
# Input Format:
#   Two-column tab-separated file:
#   <GO_term>  <p-value>
#
# Author: htafer
# Last Updated: 2025-07-28
# ====================================

use strict;
use warnings;
use autodie;
use feature 'say';

use Data::Dumper;
use WWW::Mechanize;
use WWW::Mechanize::FormFiller;
use URI::URL;
use File::Slurp;
use Getopt::Std;
use Log::Simple;
use Try::Tiny;

# Constants
use constant {
    DEFAULT_CUTOFF => 0.7,
    DEFAULT_IS_PVALUE => 'yes',
    DEFAULT_ORDER => 'higher',
    DEFAULT_UNIVERSE => 4932,
    DEFAULT_MEASURE => 'SIMREL',
    DEFAULT_TYPE => 'BP',
    REVIGO_URL => 'http://revigo.irb.hr/',
};

# Initialize logging
my $log = Log::Simple->new();

# Process command line arguments
my %opts;
getopts("f:c:i:w:g:m:t", \%opts);

# Show usage if no input file specified
if (!defined($opts{f})) {
    show_usage();
    exit 1;
}
  # Function to show usage information
sub show_usage {
    say "
Usage: revigo.pl -f <input_file> [options]

Description:
    Command-line interface to REVIGO (REduce + VIsualize Gene Ontology)
    Submits GO terms to http://revigo.irb.hr/ for redundancy reduction.

Required:
    -f <file>    Input file with GO terms and p-values (tab-separated)

Options:
    -c <float>   Clustering cutoff [default: 0.7]
    -i <yes|no>  Second column is p-value [default: yes]
    -w <string>  Order (higher/lower) [default: higher]
    -g <int>     GO universe size [default: 4932]
    -m <string>  Similarity measure [default: SIMREL]
    -t <string>  GO type (BP|CC|MF|all) [default: BP]

Example:
    revigo.pl -f go_terms.txt -c 0.8 -t BP
";
}

# Validate options and set defaults
sub validate_options {
    my ($opts) = @_;
    
    # Check input file exists and is readable
    die "Input file not found: $opts->{f}
" unless -f $opts->{f};
    die "Input file not readable: $opts->{f}
" unless -r $opts->{f};
    
    # Set defaults for optional parameters
    $opts->{c} = DEFAULT_CUTOFF unless defined $opts->{c};
    $opts->{i} = DEFAULT_IS_PVALUE unless defined $opts->{i};
    $opts->{w} = DEFAULT_ORDER unless defined $opts->{w};
    $opts->{g} = DEFAULT_UNIVERSE unless defined $opts->{g};
    $opts->{m} = DEFAULT_MEASURE unless defined $opts->{m};
    $opts->{t} = DEFAULT_TYPE unless defined $opts->{t};
    
    # Validate parameter values
    die "Invalid cutoff value: $opts->{c}
" unless $opts->{c} =~ /^\d*\.?\d+$/ && $opts->{c} > 0 && $opts->{c} <= 1;
    die "Invalid p-value flag: $opts->{i}
" unless $opts->{i} =~ /^(yes|no)$/i;
    die "Invalid order: $opts->{w}
" unless $opts->{w} =~ /^(higher|lower)$/i;
    die "Invalid GO universe size: $opts->{g}
" unless $opts->{g} =~ /^\d+$/;
    die "Invalid GO type: $opts->{t}
" unless $opts->{t} =~ /^(BP|CC|MF|all)$/i;
}

# Process GO terms through REVIGO
sub process_revigo {
    my ($opts) = @_;
    
    $log->info("Reading input file: $opts->{f}");
    my $go_string = try {
        read_file($opts->{f})
    } catch {
        die "Error reading input file: $_
"
    };
    
    $log->info("Connecting to REVIGO...");
    my $agent = WWW::Mechanize->new(
        autocheck => 1,
        timeout => 60,
        agent => 'RevigoClient/1.0'
    );
    my $formfiller = WWW::Mechanize::FormFiller->new();
    
    try {
        # Configure proxy if needed
        $agent->env_proxy();
        
        # Connect to REVIGO
        $agent->get(REVIGO_URL);
        die "Could not connect to REVIGO
" unless $agent->success;
        
        # Find and fill form
        $agent->form_number(1) if $agent->forms and scalar @{$agent->forms};
        die "Could not find REVIGO submission form
" unless $agent->current_form;
        
        # Fill form fields
        $formfiller->add_filler('goList' => Fixed => $go_string);
        $formfiller->add_filler('cutoff' => Fixed => $opts->{c});
        $formfiller->add_filler('isPValue' => Fixed => $opts->{i});
        $formfiller->add_filler('whatIsBetter' => Fixed => $opts->{w});
        $formfiller->add_filler('goSizes' => Fixed => $opts->{g});
        $formfiller->add_filler('measure' => Fixed => $opts->{m});
        $formfiller->fill_form($agent->current_form);

        # Submit form
        $log->info("Submitting data to REVIGO...");
        $agent->click("startRevigo");
        die "Form submission failed\n" unless $agent->success;
        
        # Get results based on GO type
        $log->info("Retrieving results for GO type: $opts->{t}");
        my $table_num = {
            'BP' => 1,
            'CC' => 2,
            'MF' => 3
        }->{uc($opts->{t})};
        
        if ($table_num) {
            my $link = $agent->follow_link(url => "toR_treemap.jsp?table=$table_num");
            die "Could not retrieve results\n" unless $link && $link->{_content};
            print $link->{_content};
        } elsif (lc($opts->{t}) eq 'all') {
            # Handle 'all' option by getting all tables
            for my $type (qw(BP CC MF)) {
                my $num = {
                    'BP' => 1,
                    'CC' => 2,
                    'MF' => 3
                }->{$type};
                
                $log->info("Retrieving results for $type...");
                my $link = $agent->follow_link(url => "toR_treemap.jsp?table=$num");
                die "Could not retrieve results for $type\n" unless $link && $link->{_content};
                print "# $type Results\n", $link->{_content}, "\n";
            }
        }
    } catch {
        $log->error("Error processing REVIGO request: $_");
        die "REVIGO processing failed: $_\n";
    };
}

# Main execution
try {
    # Validate options
    validate_options(\%opts);
    
    # Process through REVIGO
    process_revigo(\%opts);
    
    $log->info("Processing completed successfully");
    exit 0;
} catch {
    $log->error("Fatal error: $_");
    die "Fatal error: $_\n";
}
