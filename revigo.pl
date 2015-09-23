#!/usr/bin/perl -w

use Data::Dumper;
use strict;
use warnings;
use WWW::Mechanize;
use WWW::Mechanize::FormFiller;
use URI::URL;
use File::Slurp;
use Getopt::Std;

my %opts;
getopt("f:c:i:w:g:m:t",\%opts);

my $go_string;
if(!defined($opts{f})){
  print"
revigo.pl is a cli interface to the revigo.irb.hr website\n
-f [string]  file containing the GO p-value entries\n
-c [decimal]{0.7} cutoff on cluster score\n
-i [yes/no]{yes}  is the second entry a p-value\n
-w [string]{higher} is the order higher is better\n
-g [number]{4932} the size of the organism go universe\n
-m [string]{SIMREL} the algorithm to compute the similarity score\n
-t [string]{BP} the GO type (BP|CC|MF|all)
revigo.pl $opts{f}";
}
else{
  $go_string = read_file("$opts{f}");
  if(!defined($opts{c})){$opts{c}=0.7;}
  if(!defined($opts{i})){$opts{i}='yes';}
  if(!defined($opts{w})){$opts{w}='higher';}
  if(!defined($opts{g})){$opts{g}='4932';}
  if(!defined($opts{m})){$opts{m}='SIMREL';}
  if(!defined($opts{t})){$opts{m}='BP';}
  my $agent = WWW::Mechanize->new( autocheck => 1 );
  my $formfiller = WWW::Mechanize::FormFiller->new();
  $agent->env_proxy();
  $agent->get('http://revigo.irb.hr/');
  $agent->form_number(1) if $agent->forms and scalar @{$agent->forms};
  $formfiller->add_filler( 'goList' => Fixed => $go_string);
  $formfiller->add_filler( 'cutoff' => Fixed => $opts{c} );
  $formfiller->add_filler( 'isPValue' => Fixed => $opts{i} );
  $formfiller->add_filler( 'whatIsBetter' => Fixed => $opts{w} );
  $formfiller->add_filler( 'goSizes' => Fixed => $opts{g} );
  $formfiller->add_filler( 'measure' => Fixed => 'SIMREL' );
  $formfiller->fill_form($agent->current_form);

  $agent->click("startRevigo");
  #$agent->follow_link(url => 'revigo.jsp#fragment-1c');
  if($opts{t} eq 'BP'){
    my $link = $agent->follow_link(url => 'toR_treemap.jsp?table=1');
    print $link->{_content};
  }
  elsif($opts{t} eq 'CC'){
    my $link = $agent->follow_link(url => 'toR_treemap.jsp?table=2');
    print $link->{_content};
  }
  elsif($opts{t} eq 'MF'){
    my $link = $agent->follow_link(url => 'toR_treemap.jsp?table=3');
    print $link->{_content};
  }
}
