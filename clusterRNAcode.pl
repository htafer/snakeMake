#!/usr/bin/perl                                                                                                                                                                                                                                                                  
# cat rnaCode.xx* | sort -k 1,1 -k 4,4g -k 5,5g | grep Supercontig | perl ~/bin/clusterRNAcode.pl
use warnings;
use strict;
use Data::Dumper;

(my $prChr, my $prStart, my $prStop,my $prStrand)=("",-1,-1,"/");
my $locus;
$locus->{prChr}="";
$locus->{prStart}=-1;
$locus->{prStop}=-1;
$locus->{prStrand}="/";
$locus->{prScore}="1";
$locus->{prPValue}="1";
$locus->{count}=0;
while(my $line=<STDIN>){
    my @data = split(/\t/,$line);
    (my $pValue, my $score) = split(/\|/,$data[5]);
    if($data[0] ne $locus->{prChr} ||
       $data[3] >  $locus->{prStop} ||
       $data[6] ne $locus->{prStrand}){
        print $locus->{prChr},"\tRNAcode\tCDS\t",$locus->{prStart},"\t",$locus->{prStop},"\t",$locus->{prScore},"|",$locus->{prPValue},"\t",$locus->{prStrand},"\t.\t","gene_id \"Gene",$locus->{count},"\"\n";
        $locus->{count}++;

        $locus->{prChr}   = $data[0];
        $locus->{prStart} = $data[3];
        $locus->{prStop}  = $data[4];
        $locus->{prStrand} = $data[6];
        $locus->{prScore} = $score;
        $locus->{prPValue}= $pValue;
    }
    else{

        $locus->{prStart}=($locus->{prStart}>$data[3]?$data[3]        :$locus->{prStart});
        $locus->{prStop} =($locus->{prStop}>$data[4] ?$locus->{prStop}:$data[4]         );
        $locus->{prPValue} = ($locus->{prPValue} < $pValue ? $locus->{prPValue} : $pValue);
        $locus->{prScore} = ($locus->{prScore} < $score ? $locus->{prScore} : $score);
    }
}

print $locus->{prChr},"\tRNAcode\tCDS\t",$locus->{prStart},"\t",$locus->{prStop},"\t",$locus->{prScore},"|",$locus->{prPValue},"\t",$locus->{prStrand},"\t.\t","gene_id \"Gene",$locus->{count},"\"\n";


