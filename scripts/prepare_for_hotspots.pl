#!/usr/bin/perl

# Small script utility 
# Could be done in awk or bedtools and integrated in the main script directly

# Takes as input the sorted output of bedtools intersect between a gapless genome file and a coverage file
# Outputs the coverage on each segment, to be averaged later on by  bedtools and then a R script.

use POSIX qw/floor/;
use Getopt::Long;

GetOptions('win=i' => \$winsize);

while(<>){
  chomp;
  @F=split;
  $b=floor($F[2]/$winsize);
  if(floor($F[1]/$winsize) !=$b ){
    print "$F[0]\t$F[1]\t",$b*$winsize,"\t$F[3]\t",$b*$winsize-$F[1],"\t",($b*$winsize-$F[1])*$F[3],"\n";
    print "$F[0]\t",$b*$winsize,"\t$F[2]\t$F[3]\t",$F[2]-$b*$winsize,"\t",($F[2]-$b*$winsize)*$F[3],"\n";
  }else{
    $l=$F[2]-$F[1];
    print;
    print "\t$l\t",$l*$F[3],"\n";
  }
}




