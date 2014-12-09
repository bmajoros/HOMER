#!/usr/bin/perl
use strict;
$|=1;

my $files=`ls chunks`;
my @files=split/\s+/,$files;
my ($wins,$losses,$totalBaum,$totalNull);
foreach my $file (@files)
  {
    next unless($file=~/\.fasta/);

    my $baumScore=0+`/home/bmajoros/GF/homer/score-best-path 0th-ord-introns.hmm introns/chunks/$file`;
    my $nullScore=0+`/home/bmajoros/GF/homer/score-best-path 0th-ord-scrambled.hmm introns/chunks/$file`;

    #my $baumScore=0+`/home/bmajoros/GF/homer/score-best-path 5th-ord-labeled-chain.hmm scrambled/$file`;
    #my $nullScore=0+`/home/bmajoros/GF/homer/score-best-path 5th-ord-scrambled-chain.hmm scrambled/$file`;

    #my $baumScore=0+`/home/bmajoros/GF/homer/score-best-path 5th-ord-intron5.hmm introns/chunks/$file`;
    #my $nullScore=0+`/home/bmajoros/GF/homer/score-best-path 5th-ord-scrambled5.hmm introns/chunks/$file`;

    $totalNull+=$nullScore;
    $totalBaum+=$baumScore;
    print "$baumScore\t$nullScore\t";
    if($baumScore>$nullScore) {++$wins; print "\t\t*****"} else {++$losses}
    print "\n";
  }
my $success=int($wins/($wins+$losses)*100+0.5);
print "$wins, $success\%\n";






