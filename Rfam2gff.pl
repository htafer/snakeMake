#!/usr/bin/perl -w

use strict;
use warnings;

# perl ./trnascan2gff.pl < tRNAscanSE.out
# cat tRNAscanSE.out | perl ./trnascan2gff.pl
# chmod 755 ./trnascan2gff.pl < input
# ./trnascan2gff.pl < input

my $foundHash;

while(my $line=<STDIN>){
    next if($line=~/^#/);
    my @data = split(/\s+/,$line);
    next if($data[15]>0.001);
    if(defined($foundHash->{$data[2]})){
        $foundHash->{$data[2]}++;
    }
    else{
        $foundHash->{$data[2]}=1;
    }
    if($data[7] < $data[8]){
        print "$data[0]\tRfam\tgene\t$data[7]\t$data[8]\t$data[14]\t+\t.\tID=TU.$data[0].$data[2].$foundHash->{$data[2]};Name=TU.$data[0].$data[2].$foundHash->{$data[2]}\n";
    }
    else {
        print "$data[0]\tRfam\tgene\t$data[8]\t$data[7]\t$data[14]\t-\t.\tID=TU.$data[0].$data[2].$foundHash->{$data[2]};Name=TU.$data[0].$data[2].$foundHash->{$data[2]}\n";
    }
}
