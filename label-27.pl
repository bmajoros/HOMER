#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use FastaReader;
use FastaWriter;
$|=1;

$0=~/([^\/]+)$/;
my $program=$1;
my $usage="$program <in.gff> <in.fasta> <margin> <out.labeled>";
die "$usage\n" unless @ARGV==4;
my ($inGff,$inFasta,$margin,$out)=@ARGV;

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
	    $labels[$j]=$state;
	    $phase=($phase+1)%3;
	  }
	if($i+1<$numExons)
	  { # process intron
	    my $begin=$exon->getEnd();
	    my $end=$transcript->getIthExon($i+1)->getBegin();
	    my $delta;#=5*($phase-1)%3;              # SIZE OF INTRON MODULE
	    if($phase==1)    {$delta=0}
	    elsif($phase==2) {$delta=5}
	    else             {$delta=10}
	    for(my $j=$begin+2 ; $j<$end-2 ; ++$j)
	      {$labels[$j]=15+$delta}                          # INTRON STATE
	    $labels[$begin]=13+$delta;                         #  DONOR STATE
	    $labels[$begin+1]=14+$delta;                       #  DONOR STATE
            $labels[$end-2]=16+$delta;                       # ACCEPTOR STATE
	    $labels[$end-1]=17+$delta;                       # ACCEPTOR STATE
	  }
      }

    # Process start codon
    my $begin=$transcript->getBegin();
    my $startCodon=substr($subseq,$begin,3);
    if($startCodon ne "ATG") {die "bad start codon: $startCodon\n"}
    $labels[$begin]=2;                                          # START CODON
    $labels[$begin+1]=3;                                        # START CODON
    $labels[$begin+2]=4;                                        # START CODON

    # Process stop codon
    my $end=$transcript->getEnd()-3;
    my $stopCodon=substr($subseq,$end,3);
    $labels[$end]=8;                                                   # STOP
    if($stopCodon eq "TAG")    {$labels[$end+1]=9;  $labels[$end+2]=10}   # STOP
    elsif($stopCodon eq "TGA") {$labels[$end+1]=11; $labels[$end+2]=12}   # STOP
    elsif($stopCodon eq "TAA") {$labels[$end+1]=9;  $labels[$end+2]=12}   # STOP
    else {die "bad stop codon: $stopCodon\n"}

    # Output labeled sequence
    for(my $i=0 ; $i<$subseqLen ; ++$i)
      {
	my $letter=substr($subseq,$i,1);
	my $label=$labels[$i];
	print OUT "$letter $label\n";
      }
    print OUT "\n"; # blank line = end of sequence
  }






