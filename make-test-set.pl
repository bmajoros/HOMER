#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;
use GffTranscriptReader;
use Translation;
$|=1;

my $usage="$0 <*.gff> <*.fasta> <max#> <margin-size>";
die "$usage\n" unless @ARGV==4;
my ($gffFile,$fastaFile,$maxN,$margin)=@ARGV;

if(!-e "chunks") {system("mkdir chunks")}

##################### LOAD GFF #######################
my $gffReader=new GffTranscriptReader;
my $transcripts=$gffReader->loadGFF($gffFile);
my $numTranscripts=@$transcripts;

# Index by contig ID and add padding around each gene
#if($numTranscripts>$maxN) {$numTranscripts=$maxN}
my %records;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $begin=$transcript->getBegin();
    my $end=$transcript->getEnd();
    die unless $begin<$end;
    next unless $begin>1;
    $begin-=$margin;
    $end+=$margin;
    my $substrate=$transcript->getSubstrate();
    my $strand=$transcript->getStrand();
    push @{$records{"\L$substrate"}},[$begin,$end-$begin,$strand,$transcript];
  }

# Eliminate any gene having alternative splicing
my @substrates=keys %records;
my $n=@substrates;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $substrate=$substrates[$i];
    my @mask;
    my $recordsOnSubstrate=$records{$substrate};
    next unless defined $recordsOnSubstrate;
    my $nn=@$recordsOnSubstrate;
    for(my $j=0 ; $j<$nn ; ++$j)
      {
	my $record=$recordsOnSubstrate->[$j];
	my $transcript=$record->[3];
	my $begin=$transcript->getBegin();
	my $end=$transcript->getEnd();
	my $len=$end-$begin;
	for(my $k=$begin ; $k<$end ; ++$k)
	  {
	    if(defined $mask[$k]) 
	      {
		$transcript->{alt}=1;
		$mask[$k]->{alt}=1;
	      }
	    else {$mask[$k]=$transcript}
	  }
      }
  }

################ PROCESS FASTA FILE #################
my $fastaReader=new FastaReader($fastaFile);
my $fastaWriter=new FastaWriter;
my $nextId=1;
my %enumeration;
my %debugHash;
my $count;
while(1)
  {
    last unless $count<$maxN;
    my ($defline,$seq)=$fastaReader->nextSequenceRef();
    last unless defined $defline;
    my $seqLen=length $$seq;
    $defline=~/>(\S+)/ || die $defline;
    my $bac="\L$1";
    my $records=$records{$bac};
    next unless defined $records;
    my $n=@$records;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	last unless $count<$maxN;
	my $record=$records->[$i];
	my ($begin,$len,$strand,$transcript)=@$record;
	next if($transcript->{alt});
	#$len+=3; ### sometimes the GFF coords don't include the stop codon...

	if($begin<0) {$begin=0}
	my $end=$begin+$len;#+3;
	if($end>$seqLen) {$end=$seqLen}
	$len=$end-$begin;

	my $transcriptId=$transcript->getID();
	if($transcriptId=~/(.*);/) {$transcriptId=$1}
	if(!defined($enumeration{$transcriptId})) 
	  {$enumeration{$transcriptId}=$nextId++}

	$debugHash{$enumeration{$transcriptId}}=$transcriptId;
	my $d=$enumeration{$transcriptId};

	$transcriptId=$enumeration{$transcriptId};
	my $outfasta="$transcriptId.fasta";
	my $outgff="$transcriptId.gff";
	open(OUTFASTA,">chunks/$outfasta") || die "can't create $outfasta";
	open(OUTGFF,">chunks/$outgff") || die "can't create $outgff";

	#print "b=$begin len=$len\n";
	my $string=substr($$seq,$begin,$len);
	my $substrate=$transcript->getSubstrate();
	$fastaWriter->addToFasta(">$transcriptId ($substrate)",$string,
				 \*OUTFASTA);
	$transcript->setSubstrate($transcriptId);

	### sometimes the GFF coords don't include the stop codon...
	my $numExons=$transcript->numExons();
	my $lastExon=$transcript->getIthExon($numExons-1);
	my $lastExonEnd=$lastExon->getEnd();
	if($strand eq "+")
	  {
	    my $stopCodonBegin=$lastExonEnd-3;
	    my $stopCodon=substr($$seq,$stopCodonBegin,3);
	    if($lastExon->getType eq "final-exon" &&
	       $stopCodon ne "TAG" && $stopCodon ne "TAA" 
	       && $stopCodon ne "TGA")
	      {
		$stopCodonBegin+=3;
		$lastExon->setEnd($lastExonEnd+3);
		$stopCodon=substr($$seq,$stopCodonBegin,3);
		if($stopCodon ne "TAG" && $stopCodon ne "TAA" && 
		   $stopCodon ne "TGA")
		  {
		    print "WARNING!  Transcript \"$transcriptId\" " .
		      "on $strand strand has no stop codon\n";
		    $lastExon->setEnd($lastExonEnd);
		   }
	      }
	  }
	else # $strand eq "-"
	  {
	    my $stopCodonBegin=$lastExon->getBegin();
	    my $stopCodon=substr($$seq,$stopCodonBegin,3);
	    $stopCodon=Translation::reverseComplement(\$stopCodon);
	    if($lastExon->getType() eq "final-exon" &&
	       $stopCodon ne "TAG" && $stopCodon ne "TAA" 
	       && $stopCodon ne "TGA")
	      {
		$stopCodonBegin-=3;
		$lastExon->setBegin($stopCodonBegin);
		$stopCodon=substr($$seq,$stopCodonBegin,3);
		$stopCodon=Translation::reverseComplement(\$stopCodon);
		if($stopCodon ne "TAG" && $stopCodon ne "TAA" && 
		   $stopCodon ne "TGA")
		  {
		    my $d=$debugHash{$transcriptId};
		    print "WARNING!  Transcript \"$transcriptId\" " .
		      "($d) on $strand strand has no stop codon\n";
		    $lastExon->setBegin($stopCodonBegin+3);
		  }
	      }
	  }
	
	# Shift coords and output GFF
	my $delta=-$begin;
	for(my $j=0 ; $j<$numExons ; ++$j)
	  {
	    my $exon=$transcript->getIthExon($j);
	    $exon->shiftCoords($delta);
	    my $gff=$exon->toGff();
	    print OUTGFF "$gff";
	  }
	close(OUTGFF);
	close(OUTFASTA);
	++$count;
      }
  }

