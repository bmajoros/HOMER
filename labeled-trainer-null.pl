#!/usr/bin/perl
use strict;
use NgramIterator;
$|=1;

$0=~/([^\/]+)$/;
my $usage="$1 <labeled-sequence> <highest-state-num> <order> <out.hmm>";
die "$usage\n" unless @ARGV==4;
my ($seqFile,$highestState,$order,$outfile)=@ARGV;
my $MIN_SAMPLE_SIZE=175;

# Initialization
my (%symbolMap,$label,$symbolId,@baseAlphabet);
@baseAlphabet=("A","C","G","N","T");
$symbolId=0;
for(my $i=1 ; $i<=$order+1 ; ++$i)
  {
    my $iterator=new NgramIterator("ACGNT",$i);
    while($label=$iterator->nextString())
      {
	$symbolMap{$label}=$symbolId;
	++$symbolId;
      }
  }
my $NUM_SYMBOLS=$symbolId;

# Read input and gather statistics
my (@transitions,@emissions,$prevState);
open(IN,$seqFile) || die "Can't open file: $seqFile\n";
while(<IN>)
  {
    if(/^([ACGTN]+)\s+(\d+)/)
      {
	my ($ngram,$state)=($1,$2);
	my $len=length $ngram;
	for(my $i=0 ; $i<$len ; ++$i)
	  {++$emissions[$state][$symbolMap{substr($ngram,$i,$len)}]}
	if(defined $prevState) {++$transitions[$prevState][$state]}
	else {++$transitions[0][$state]}
	$prevState=$state;
      }
    elsif(/^\s*$/)
      {
	++$transitions[$prevState][0];
	undef $prevState;
      }
  }
close(IN);

# Convert integral counts into probabilities
for(my $state=0 ; $state<=$highestState ; ++$state)
  {
    my $totalTrans=0;
    for(my $nextState=0 ; $nextState<=$highestState ; ++$nextState)
      {	$totalTrans+=$transitions[$state][$nextState] }
    if($totalTrans)
      {	for(my $nextState=0 ; $nextState<=$highestState ; ++$nextState)
	  { $transitions[$state][$nextState]/=$totalTrans }
      }
    next unless $state;

    for(my $i=0 ; $i<=$order ; ++$i)
      {
	my $iterator=new NgramIterator("ACGNT",$i);
	while(1)
	  {
	    my $history=$iterator->nextString();
	    last unless defined $history;
	    my $sum=0;
	    foreach my $base (@baseAlphabet)
	      {
		my $label="$history$base";
		$sum+=$emissions[$state][$symbolMap{$label}];
	      }
	    foreach my $base (@baseAlphabet)
	      {
		my $label="$history$base";
		my $symbolId=$symbolMap{$label};
		if($sum==0 && $i==0) {die "0 sum for state $state"}
		if($sum>=$MIN_SAMPLE_SIZE || $i==0)
		  {$emissions[$state][$symbolId]/=$sum}
		else
		  {undef $emissions[$state][$symbolId]}
	      }
	  }
      }
  }

# Output HMM
open(OUT,">$outfile") || die "Can't write to file: $outfile\n";
my $numStates=$highestState+1;
print OUT "$numStates $order\n";
for(my $i=0 ; $i<$numStates ; ++$i)
  {
    for(my $j=0 ; $j<$numStates ; ++$j)
      {
	my $p=int(10000*$transitions[$j][$i]+0.5)/10000;
	print OUT "$p\t";
      }
    print OUT "\n";
  }
for(my $i=0 ; $i<$NUM_SYMBOLS ; ++$i)
  {
    for(my $j=0 ; $j<$numStates ; ++$j)
      {
	my $p=int(10000*$emissions[$j][$i]+0.5)/10000;
	print OUT "$p\t";
      }
    print OUT "\n";
  }

close(OUT);


