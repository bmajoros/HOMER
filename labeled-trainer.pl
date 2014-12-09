#!/usr/bin/perl
use strict;

$0=~/([^\/]+)$/;
my $usage="$1 <labeled-sequence> <highest-state-num> <out.hmm>";
die "$usage\n" unless @ARGV==3;
my ($seqFile,$highestState,$outfile)=@ARGV;

# Initialization
my %symbolMap;
my @symbols=('A','C','G','N','T');
for(my $i=0 ; $i<@symbols ; ++$i) {$symbolMap{$symbols[$i]}=$i}

# Read input and gather statistics
my (@transitions,@emissions,$prevState);
open(IN,$seqFile) || die "Can't open file: $seqFile\n";
while(<IN>)
  {
    if(/^([ACGTN])\s+(\d+)/)
      {
	my ($symbol,$state)=($1,$2);
	++$emissions[$state][$symbolMap{$symbol}];
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

    my $totalEmit=0;
    for(my $i=0 ; $i<5 ; ++$i)
      {	$totalEmit+=$emissions[$state][$i] }
    if($totalEmit)
      {	for(my $i=0 ; $i<5 ; ++$i)
	  { $emissions[$state][$i]/=$totalEmit }
      }
  }

# Output HMM
open(OUT,">$outfile") || die "Can't write to file: $outfile\n";
my $numStates=$highestState+1;
print OUT "$numStates\n";
for(my $i=0 ; $i<$numStates ; ++$i)
  {
    for(my $j=0 ; $j<$numStates ; ++$j)
      {
	my $p=0+$transitions[$j][$i];
	print OUT "$p\t";
      }
    print OUT "\n";
  }
for(my $i=0 ; $i<5 ; ++$i)
  {
    for(my $j=0 ; $j<$numStates ; ++$j)
      {
	my $p=0+$emissions[$j][$i];
	print OUT "$p\t";
      }
    print OUT "\n";
  }

close(OUT);


