#!/usr/bin/perl -w

use strict;
use warnings;

# perl ./trnascan2gff.pl < tRNAscanSE.out
# cat tRNAscanSE.out | perl ./trnascan2gff.pl
# chmod 755 ./trnascan2gff.pl < input
# ./trnascan2gff.pl < input

<STDIN>;
<STDIN>;
<STDIN>;

while(my $line=<STDIN>){
    my @data = split(/\s+/,$line);
    if($data[2] < $data[3]){
	#print gene
        print "$data[0]\ttRNAscan-SE\tgene\t$data[2]\t$data[3]\t$data[8]\t+\t.\tID=TU.$data[0].tRNA-$data[4].$data[1];Name=TU.$data[0].tRNA-$data[4].$data[1]\n";
        #print exon
        if($data[6]+$data[7]>0){
            print "$data[0]\ttRNAscan-SE\texon\t$data[2]\t$data[6]\t$data[8]\t+\t.\tID=exon.$data[0].tRNA-$data[4].$data[1].exon1; Parent=TU.$data[0].tRNA-$data[4].$data[1]\n";
            print "$data[0]\ttRNAscan-SE\texon\t$data[7]\t$data[3]\t$data[8]\t+\t.\tID=exon.$data[0].tRNA-$data[4].$data[1].exon2; Parent=TU.$data[0].tRNA-$data[4].$data[1]\n";
        }
        else{
            print "$data[0]\ttRNAscan-SE\texon\t$data[2]\t$data[3]\t$data[8]\t+\t.\tID=exon.$data[0].tRNA-$data[4].$data[1].exon1; Parent=TU.$data[0].tRNA-$data[4].$data[1]\n";
        }
    }
    else {
        #print gene
        print "$data[0]\ttRNAscan-SE\tgene\t$data[3]\t$data[2]\t$data[8]\t-\t.\tID=TU.$data[0].tRNA-$data[4].$data[1];Name=TU.$data[0].tRNA-$data[4].$data[1]\n";
        #print exon
        if($data[6]+$data[7]>0){
            print "$data[0]\ttRNAscan-SE\texon\t$data[3]\t$data[7]\t$data[8]\t-\t.\tID=exon.$data[0].tRNA-$data[4].$data[1].exon1; Parent=TU.$data[0].tRNA-$data[4].$data[1]\n";
            print "$data[0]\ttRNAscan-SE\texon\t$data[6]\t$data[2]\t$data[8]\t-\t.\tID=exon.$data[0].tRNA-$data[4].$data[1].exon2; Parent=TU.$data[0].tRNA-$data[4].$data[1]\n";
        }
        else{
            print "$data[0]\ttRNAscan-SE\texon\t$data[3]\t$data[2]\t$data[8]\t-\t.\tID=exon.$data[0].tRNA-$data[4].$data[1].exon1; Parent=TU.$data[0].tRNA-$data[4].$data[1]\n";
        }
    }
}
