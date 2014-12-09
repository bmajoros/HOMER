#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use FastaReader;
use FastaWriter;
$|=1;

$0=~/([^\/]+)$/;
my $program=$1;
my $usage="$program <in.gff> <in.fasta> <margin> <order> <out.labeled>";
die "$usage\n" unless @ARGV==5;
my ($inGff,$inFasta,$margin,$order,$out)=@ARGV;

my $gffReader=new GffTranscriptReader();
my $genesByContig=$gffReader->hashBySubstrate($inGff);

open(OUT,">$out") || die "Can't create file: $out\n";
my $fastaReader=new FastaReader($inFasta);
my $fastaWriter=new FastaWriter;
while(1)
  {
    my ($def,$seq)=$fastaReader->nextSequence();
    last unless defined $def;
    $def=~/^\s*>\s*(\S+)/ || die "can't parse defline: $def\n";
    my $contigId=$1;
    my $genesOnContig=$genesByContig->{$contigId};
    next unless defined $genesOnContig;
    my $numGenes=@$genesOnContig;
    for(my $i=0 ; $i<$numGenes ; ++$i)
      {
	my $gene=$genesOnContig->[$i];
	process($gene,\$seq);
      }
  }
close(OUT);

#=============================================================
sub process
  {
    my ($transcript,$seqRef)=@_;

    # Extract subsequence to process
    my $seqLength=length $$seqRef;
    my $begin=$transcript->getBegin()-$margin;
    if($begin<0) {$begin=0}
    my $end=$transcript->getEnd()+$margin;
    if($end>=$seqLength) {$end=$seqLength}
    my $subseqLen=$end-$begin;
    my $subseq=substr($$seqRef,$begin,$subseqLen);
    $transcript->shiftCoords(-$begin);
    my @labels;

    # First, label everything intergenic as default
    for(my $i=0 ; $i<$subseqLen ; ++$i)
      {$labels[$i]=1}                                      # INTERGENIC STATE

    # Process all exons and introns
    my $phase=0;
    my $numExons=$transcript->numExons();
    for(my $i=0 ; $i<$numExons ; ++$i)
      {
	my $exon=$transcript->getIthExon($i);
	my $begin=$exon->getBegin();
	my $end=$exon->getEnd();
	for(my $j=$begin ; $j<$end ; ++$j)
	  { # process exon
	    my $state=5+$phase;                              # EXON STATES
	    $labels[$j]=$state
	      unless($j-$begin<3 && $i>0);
	    $phase=($phase+1)%3;
	  }
	if($i+1<$numExons)
	  { # process intron
	    my $begin=$exon->getEnd();
	    my $end=$transcript->getIthExon($i+1)->getBegin();
	    my $delta;#=5*($phase-1)%3;            # SIZE OF INTRON MODULE
	    if($phase==1)    {$delta=0}
	    elsif($phase==2) {$delta=5}
	    else             {$delta=10}
	    for(my $j=$begin+2 ; $j<$end-2 ; ++$j)
	      {$labels[$j]=15+$delta}                          # INTRON STATE
	    $labels[$begin]=13+$delta;                         #  DONOR STATE
	    $labels[$begin+1]=14+$delta;                       #  DONOR STATE
            $labels[$end-2]=16+$delta;                       # ACCEPTOR STATE
	    $labels[$end-1]=17+$delta;                       # ACCEPTOR STATE

	    # Process 3 bases before donor
	    my $d;
	    if($phase==1)    {$d=0}
	    elsif($phase==2) {$d=6}
	    else             {$d=12}
	    $labels[$begin-3]=78+$d;
	    $labels[$begin-2]=79+$d;
	    $labels[$begin-1]=80+$d;

	    # Process 3 bases after acceptor
	    $labels[$end]=81+$d;
	    $labels[$end+1]=82+$d;
	    $labels[$end+2]=83+$d;

	    # Process donor & acceptor context windows
	    $delta*=2;
	    if($end-$begin>15)
	      { for(my $i=0 ; $i<5 ; ++$i)
		  { $labels[$begin+2+$i]=48+$delta+$i} # DONOR CONTEXT WINDOW
		for(my $i=1 ; $i<=5 ; ++$i)
		  { $labels[$end-2-$i]=58+$delta-$i}# ACCEPTOR CONTEXT WINDOW
	      }
	  }
      }

    # Process start codon
    my $begin=$transcript->getBegin();
    my $startCodon=substr($subseq,$begin,3);
    if($startCodon ne "ATG") {die "bad start codon: $startCodon\n"}
    $labels[$begin]=2;                                          # START CODON
    $labels[$begin+1]=3;                                        # START CODON
    $labels[$begin+2]=4;                                        # START CODON

    # Process upstream UTR
    my $utrEnd=$begin;
    my $utrBegin=$utrEnd-10;
    if($utrBegin>=0)
      {	for(my $i=$utrBegin ; $i<$utrEnd ; ++$i)
	  { $labels[$i]=28+$i-$utrBegin }               # UPSTREAM UTR STATES
      }
    $utrEnd=$utrBegin;
    $utrBegin-=50;
    if($utrBegin<50) {$utrBegin=50}
    #for(my $i=$utrBegin ; $i<$utrEnd ; ++$i)
      #{ $labels[$i]=28 }                                 # UPSTREAM UTR STATE

    # Process stop codon
    my $end=$transcript->getEnd()-3;
    my $stopCodon=substr($subseq,$end,3);
    $labels[$end]=8;                                                   # STOP
    if($stopCodon eq "TAG")    {$labels[$end+1]=9;  $labels[$end+2]=10}   # STOP
    elsif($stopCodon eq "TGA") {$labels[$end+1]=11; $labels[$end+2]=12}   # STOP
    elsif($stopCodon eq "TAA") {$labels[$end+1]=9;  $labels[$end+2]=12}   # STOP
    else {die "bad stop codon: $stopCodon\n"}

    # Process downstream UTR
    my $utrBegin=$transcript->getEnd();
    my $utrEnd=$utrBegin+10;
    if($utrEnd<=$subseqLen)
      {	for(my $i=$utrBegin ; $i<$utrEnd ; ++$i)
	  { $labels[$i]=38+$i-$utrBegin }             # DOWNSTREAM UTR STATES
      }
    $utrBegin=$utrEnd;
    $utrEnd+=50;
    if($utrEnd>$subseqLen-50) {$utrEnd=$subseqLen-50}
    #for(my $i=$utrBegin ; $i<$utrEnd ; ++$i)
      #{ $labels[$i]=47 }                               # DOWNSTREAM UTR STATE

    # Output labeled sequence
    for(my $i=0 ; $i<$subseqLen ; ++$i)
      {
	#my $letter=substr($subseq,$i,1);
	my $begin=$i-$order;
	if($begin<0) {$begin=0}
	my $letter=substr($subseq,$begin,$i-$begin+1);
	my $label=$labels[$i];
	print OUT "$letter $label\n";
      }
    print OUT "\n"; # blank line = end of sequence
  }






