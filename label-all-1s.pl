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
    for(my $i=0 ; $i<$len ; ++$i)
      {
	my $begin=$i-5;
	if($begin<0) {$begin=0}
	my $ngramLen=$i-$begin+1;
	my $ngram=substr($seq,$begin,$ngramLen);
	my $state=1;
	print "$ngram $state\n";
      }
    print "\n";
  }


