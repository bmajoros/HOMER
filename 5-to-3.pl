#!/usr/bin/perl
use strict;

$0=~/([^\/]+)$/;
die "$1 <labeled-in> <labeled-out>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(IN,$infile) || die "can't open $infile\n";
open(OUT,">$outfile") || die "can't write to file: $outfile\n";
while(<IN>)
  {
    if(/(\S)\s*(\d+)/)
      {
	my ($symbol,$state)=($1,$2);
	if($state>=2 && $state<=4) {$state=2}
	elsif($state==5) {$state=3}
	print OUT "$symbol $state\n";
      }
    else {print OUT}
  }
close(OUT);
close(IN);




