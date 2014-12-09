#!/usr/bin/perl
use strict;
use FastaReader;
$|=1;

$0=~/([^\/]+)$/;
my $program=$1;
my $usage="$program <in.fasta> <order>";
die "$usage\n" unless @ARGV==2;
my ($inFasta,$order)=@ARGV;

my $fastaReader=new FastaReader($inFasta);
while(1)
  {
    my ($def,$seq)=$fastaReader->nextSequence();
    last unless defined $def;
    my $len=length $seq;
    my $skip=(rand>0.5 ? 0 : (rand>0.5 ? 1 : 2));
    for(my $i=0 ; $i<$skip ; ++$i)
      {
	my $ngramLen=$i+1;
	my $ngram=substr($seq,0,$ngramLen);
	my $state=$skip-$i;
	print "$ngram $state\n";
      }
    for(my $i=$skip ; $i<$len ; ++$i)
      {
	my $begin=$i-5;
	if($begin<0) {$begin=0}
	my $ngramLen=$i-$begin+1;
	my $ngram=substr($seq,$begin,$ngramLen);
	my $state=($i-$skip)%6+3;
	if($i==$len-2 && $state==6) {$state=9}
	elsif($i==$len-1 && $state==7) {$state=10}
	elsif($i==$len-1 && $state==6) {$state=10}
	elsif($i==$len-2 && $state==3) {$state=9}
	elsif($i==$len-1 && $state==4) {$state=10}
	elsif($i==$len-1 && $state==3) {$state=10}
	print "$ngram $state\n";
      }
    print "\n";
  }


