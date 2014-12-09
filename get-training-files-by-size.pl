#!/usr/bin/perl
use strict;
use FastaReader;
use Translation;
use GffTranscriptReader;
use CodonIterator;
use Codon;
use Exon;
use Transcript;
use FileHandle;
use FastaWriter;
$|=1;

my $INTERGENIC_LENGTH=1000;

    ##################################################################
    #
    # Note that we chop off the last two bases of any exon which is
    # followed by an intron, and the first base from any exon which
    # is preceded by an intron, because these are part of the splice-
    # junction model.
    #
    # Donor sites:     __GT___ (includes last two bases of previous exon)
    # Acceptor sites:  ____AG_ (includes first base of next exon)
    #
    #####################################################################

my $usage="$0 <*.gff> <*.fasta> <max-example-nucleotides>";
die "$usage\n" unless @ARGV==3;
my ($gffFilename,$fastaFilename,$N)=@ARGV;
$fastaFilename=~/([^\/]+)\.fasta/;
my $axisName=$1;

# Read feature coordinates from GFF file
print STDERR "Reading features from GFF file...\n";
my $gffReader=new GffTranscriptReader;
my $transcripts=$gffReader->loadGFF($gffFilename);
my $numTranscripts=@$transcripts;
my %transcripts;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $substrate=$transcript->getSubstrate();
    push @{$transcripts{$substrate}},$transcript;
  }

# Load sequence of genomic axis
print STDERR "Reading genomic axis sequence...\n";
my $fastaReader=new FastaReader($fastaFilename);

# Create output files
my $fhIntrons=new FileHandle(">introns.fasta") 
  || die "can't create introns.fasta";
my $fhExons=new FileHandle(">exons.fasta") 
  || die "can't create exons.fasta";
my $fhDonors=new FileHandle(">donors.fasta") 
  || die "can't create donors.fasta";
my $fhAcceptors=new FileHandle(">acceptors.fasta") 
  || die "can't create acceptors.fasta";
my $fhUpStreamIntergenic=new FileHandle(">upstream-intergenic.fasta") 
  || die "can't create upstream-intergenic.fasta";
my $fhDownStreamIntergenic=new FileHandle(">downstream-intergenic.fasta") 
  || die "can't create downstream-intergenic.fasta";
my $fhStartCodons=new FileHandle(">start-codons.fasta") 
  || die "can't create start-codons.fasta";
my $fhStopCodons=new FileHandle(">stop-codons.fasta") 
  || die "can't create stop-codons.fasta";
my $fhDonorFrameShifts=new FileHandle(">donor-frame-shifts.fasta") || 
  die "can't create donor-frame-shifts.fasta";
my $fhAccFrameShifts=new FileHandle(">acceptor-frame-shifts.fasta") || 
  die "can't create acceptor-frame-shifts.fasta";
my $fastaWriter=new FastaWriter;

my ($intronCounter, $exonCounter, $donorCounter, $acceptorCounter,
    $upStreamCounter, $downStreamCounter, $startCounter, $stopCounter,
    $dfsCounter, $afsCounter);

while(1)
  {
    last unless $intronCounter<$N || $exonCounter<$N || $donorCounter<$N ||
      $acceptorCounter<$N || $upStreamCounter<$N || $downStreamCounter<$N ||
      $startCounter<$N || $stopCounter<$N || $dfsCounter<$N || $afsCounter<$N;

    my ($defline,$axisSequence)=$fastaReader->nextSequence();
    last unless $defline && $axisSequence;
    $defline=~/>\s*(\S+)/ || die $defline;
    my $substrate=$1;

    my $transcripts=$transcripts{$substrate};
    next unless defined($transcripts);
    my $numTranscripts=@$transcripts;
    for(my $i=0 ; $i<$numTranscripts ; ++$i)
      {
	my $transcript=$transcripts->[$i];
	my $codonIterator=new CodonIterator($transcript,\$axisSequence);
	my $codons=$codonIterator->getAllCodons();
	my $numCodons=@$codons;
	my $startCodon=$codons->[0];
	my $stopCodon=$codons->[$numCodons-1];
	my $exons=$transcript->{exons};
	my $transcriptId=$transcript->{transcriptId};
	my $strand=$transcript->{strand};
	print STDERR "processing transcript $transcriptId on $substrate\n";
	
	# Attach codons to their exons
	attachCodonsToExons($codons);
	
	# Drop exons before start codon or after stop codon.  Also,
	# trim portion of exon before start codon, similarly for stop.
	trimExons($exons,$startCodon,$stopCodon,$transcript,$strand);
	my $numExons=@$exons;
	if($strand eq "+") 
	  {
	    $transcript->{fivePrime}=$exons->[0]->{begin};
	    $transcript->{threePrime}=$exons->[$numExons-1]->{end};
	  }
	else
	  {
	    $transcript->{fivePrime}=$exons->[$numExons-1]->{begin};
	    $transcript->{threePrime}=$exons->[0]->{end};
	  }
	
	# Process start & stop codons
	my $startCodonCoord=$startCodon->{absoluteCoord};
	if($startCounter<$N)
	  {
	    $fastaWriter->addToFasta(">$transcriptId /class=StartCodon ".
				     "/loc=$startCodonCoord /axis=$axisName",
				     $startCodon->{triplet},
				     $fhStartCodons);
	    $startCounter+=length $startCodon->{triplet};
	  }
	my $stopCodonCoord=$stopCodon->{absoluteCoord};
	if($stopCounter<$N)
	  {
	    $fastaWriter->addToFasta(">$transcriptId /class=StartCodon ".
				     "/loc=$stopCodonCoord /axis=$axisName",
				     $stopCodon->{triplet},
				     $fhStopCodons);
	    $stopCounter+=length $stopCodon->{triplet};
	  }
	
	#####################################################################
	# Process each exon and its surrounding subfeatures (introns, splice,
	# junctions, and frame shifts)
	#####################################################################
	for(my $i=0 ; $i<$numExons ; ++$i)
	  {
	    my $exon=$exons->[$i];
	    my $codons=$exon->{codons};
	    next unless defined $codons;
	    my $numCodons=@$codons;
	    
	    #####################################################
	    # Process preceding intron, acceptor, and frame shift
	    #####################################################
	    if($i>0) 
	      {
		# First, the intron:
		my $prevExon=$exons->[$i-1];
		my ($intronBegin,$intronEnd);
		if($strand eq "+")
		  {
		    $intronBegin=$prevExon->{end}+5; # skip over donor
		    $intronEnd=$exon->{begin}-6; # skip over acceptor
		  }
		else
		  {
		    $intronBegin=$exon->{end}+5; # skip over donor
		    $intronEnd=$prevExon->{begin}-6; # skip over acceptor
		  }
		my $intronLength=$intronEnd-$intronBegin;
		my $intronSeq=substr($axisSequence,$intronBegin,$intronLength);
		if($strand eq "-") 
		  {$intronSeq=Translation::reverseComplement(\$intronSeq)}
		if($intronCounter<$N)
		  {
		    $fastaWriter->addToFasta(">$transcriptId /class=intron ".
					     "/begin=$intronBegin ".
					     "/end=$intronEnd".
					     "  /axis=$axisName",
					     $intronSeq,
					     $fhIntrons);
		    $intronCounter+=length $intronSeq;
		  }
		
		# Write out acceptor site
		my $acceptorBegin=
		  ($strand eq "+" ? $exon->{begin}-6 : $exon->{end}-1);
		my $acceptorSeq=substr($axisSequence,$acceptorBegin,7);
		if($strand eq "-")
		  {$acceptorSeq=Translation::reverseComplement(\$acceptorSeq)}
		if($acceptorCounter<$N && length($acceptorSeq)>6)
		  {
		    $fastaWriter->addToFasta(">$transcriptId /class=acceptor ".
					     "/begin=$acceptorBegin ".
					     "/strand=$strand",
					     $acceptorSeq,
					     $fhAcceptors);
		    $acceptorCounter+=length $acceptorSeq;
		  }
		
		# Finally, the frame shift after the acceptor site
		my $firstCodon=$codons->[0];
		if($firstCodon->{relativeCoord}==0)
		  {
		    shift @$codons;
		    $firstCodon=$codons->[0];
		    --$numCodons;
		  }
		if(defined $firstCodon)
		  {
		    my $frameShiftBegin=
		      ($strand eq "+" ? 
		       $exon->{begin}+1 : 
		       $firstCodon->{absoluteCoord});
		    my $frameShiftEnd=
		      ($strand eq "+" ?
		       $firstCodon->{absoluteCoord} :
		       $exon->{end}-1);
		    my $frameShiftLen=$frameShiftEnd-$frameShiftBegin;
		    my $frameShiftSeq=
		      substr($axisSequence,$frameShiftBegin,$frameShiftLen);
		    if($strand eq "-")
		      {$frameShiftSeq=
			 Translation::reverseComplement(\$frameShiftSeq)}
		    if($afsCounter<$N)
		      {
			$fastaWriter->addToFasta(">$transcriptId ".
						 "/class=AcceptorFrameShift ".
						 "/begin=$frameShiftBegin ".
						 "/len=$frameShiftLen ".
						 "/strand=$strand",
						 $frameShiftSeq,
						 $fhAccFrameShifts);
			$afsCounter+=length $frameShiftSeq;
		      }
		  }
	      }
	    
	    ##################################################################
	    # Process the in-frame portion of the exon
	    ##################################################################
	    my $exonSeq;
	    my $exonLen=$exon->getLength();
	    my $firstCodonBegin;
	    for(my $i=0 ; $i<$numCodons ; ++$i)
	      {
		my $codon=$codons->[$i];
		my $relCoord=$codon->{relativeCoord};
		next if $exon->{order}>0 && $relCoord==0 || 
		  $exon->{order}<$numExons-1 && $relCoord+2>=$exonLen-2;
		$exonSeq.=$codon->{triplet};
		$firstCodonBegin=$codon->{absoluteCoord}
		  unless defined $firstCodonBegin;
	      }
	    my $codingLen=length($exonSeq);
	    if($exonCounter<$N && $codingLen>2 && $exonSeq=~/^[ATCG]+$/)
	      {
		$fastaWriter->addToFasta(">$transcriptId /class=exon ".
					 "/begin=$firstCodonBegin ".
					 "/len=$codingLen",
					 $exonSeq,
					 $fhExons);
		$exonCounter+=length $exonSeq;
	      }
	    
	    ##################################################################
	    # Process following frame shift & donor
	    ##################################################################
	    if($i<$numExons-1)
	      {
		# Write out the frame shift before the donor
		my $lastCodon=$codons->[$numCodons-1];
		die unless defined $lastCodon;
		my $exonLen=$exon->getLength();
		while($numCodons>0 && 
		      $lastCodon->{relativeCoord}+2>=$exonLen-2)
		  {
		    pop @$codons;
		    --$numCodons;
		    $lastCodon=($numCodons>=1 ? $codons->[$numCodons-1] : 
				undef);
		  }
		if(defined $lastCodon)
		  {
		    my $frameShiftBegin=
		      ($strand eq "+" ?
		       $lastCodon->{absoluteCoord}+3:
		       $exon->{begin}+2);
		    my $frameShiftEnd=
		      ($strand eq "+" ?
		       $exon->{end}-2 :
		       $lastCodon->{absoluteCoord}-3);
		    my $frameShiftLen=$frameShiftEnd-$frameShiftBegin;
		    my $frameShiftSeq=
		      substr($axisSequence,$frameShiftBegin,$frameShiftLen);
		    if($strand eq "-")
		      {$frameShiftSeq=
			 Translation::reverseComplement(\$frameShiftSeq)}
		    if($dfsCounter<$N)
		      {
			$fastaWriter->addToFasta(">$transcriptId ".
						 "/class=DonorFrameShift ".
						 "/begin=$frameShiftBegin ".
						 "/len=$frameShiftLen ".
						 "/strand=$strand",
						 $frameShiftSeq,
						 $fhDonorFrameShifts);
			$dfsCounter+=length $frameShiftSeq;
		      }
		  }
		
		# Write out donor site
		my $donorBegin=
		  ($strand eq "+" ? $exon->{end}-2 : $exon->{begin}-5);
		my $donorSeq=substr($axisSequence,$donorBegin,7);
		if($strand eq "-")
		  {$donorSeq=Translation::reverseComplement(\$donorSeq)}
		if($donorCounter<$N && length($donorSeq)>6)
		  {
		    $fastaWriter->addToFasta(">$transcriptId /class=donor ".
					     "/begin=$donorBegin ".
					     "/strand=$strand",
					     $donorSeq,
					     $fhDonors);
		    $donorCounter+=length $donorSeq;
		  }
	      }
	    undef $exon->{codons};
	  }
	
	##################################################################
	# Write out the preceding & following intergenic regions (taking
	# only the closest portions thereof)
	##################################################################
	if($i>0)
	  {
	    my $prevTranscript=$transcripts->[$i-1];
	    my ($intergenicSeq,$intergenicBegin,$intergenicEnd,$intergenicLen);
	    if($strand eq "+")
	      {
		$intergenicBegin=
		  $prevTranscript->{threePrime}+3; # skip stop codon
		$intergenicEnd=$startCodon->{absoluteCoord};
		$intergenicLen=$intergenicEnd-$intergenicBegin;
		if($intergenicLen>0)
		  {
		    $intergenicSeq=substr($axisSequence,$intergenicBegin,
					  $intergenicLen);
		  }
	      }
	    else # $strand eq "-"
	      {
		$intergenicBegin=
		  $prevTranscript->{threePrime}+3; # skip start codon
		$intergenicEnd=$transcript->{fivePrime}-3; # skip stop codon
		$intergenicLen=$intergenicEnd-$intergenicBegin;
		if($intergenicLen>0)
		  {
		    $intergenicSeq=substr($axisSequence,$intergenicBegin,
					  $intergenicLen);
		    $intergenicSeq=
		      Translation::reverseComplement(\$intergenicSeq);
		  }
	      }
	    if($intergenicLen>0)
	      {
		my $halfLength=int($intergenicLen/2);
		if($halfLength>$INTERGENIC_LENGTH) 
		  {$halfLength=$INTERGENIC_LENGTH}
		my $downStream=substr($intergenicSeq,0,$halfLength);
		my $upstreamBegin=$intergenicLen-$halfLength;
		my $upStream=substr($intergenicSeq,$upstreamBegin,$halfLength);
		my $absoluteUpstream=$intergenicEnd-$halfLength;
		if($upStreamCounter<$N && length($upStream)>=10)
		  {
		    $fastaWriter->addToFasta(">$transcriptId ".
					     "/class=upStreamIntergenic ".
					     "/begin=$absoluteUpstream ".
					     "/len=$halfLength ".
					     "/strand=$strand ",
					     $upStream,
					     $fhUpStreamIntergenic);
		    $upStreamCounter+=length $upStream;
		  }
		if($downStreamCounter<$N && length($downStream)>=10)
		  {
		    $fastaWriter->addToFasta(">$transcriptId ".
					     "/class=downStreamIntergenic ".
					     "/begin=$intergenicBegin ".
					     "/len=$halfLength ".
					     "/strand=$strand ",
					     $downStream,
					     $fhDownStreamIntergenic);
		    $downStreamCounter+=length $downStream;
		  }
	      }
	  }
	else # only one transcript on this substrate
	  {
	    my ($intergenicSeq,$intergenicBegin,$intergenicEnd,$intergenicLen);
	    if($strand eq "+")
	      {
		$intergenicBegin=0;
		$intergenicEnd=$startCodon->{absoluteCoord};
		$intergenicLen=$intergenicEnd-$intergenicBegin;
		if($intergenicLen>0)
		  {
		    $intergenicSeq=substr($axisSequence,$intergenicBegin,
					  $intergenicLen);
		  }
	      }
	    else # $strand eq "-"
	      {
		$intergenicBegin=0;
		$intergenicEnd=$transcript->{fivePrime}-3; # skip stop codon
		$intergenicLen=$intergenicEnd-$intergenicBegin;
		if($intergenicLen>0)
		  {
		    $intergenicSeq=substr($axisSequence,$intergenicBegin,
					  $intergenicLen);
		    $intergenicSeq=
		      Translation::reverseComplement(\$intergenicSeq);
		  }
	      }
	    if($intergenicLen>0)
	      {
		my $halfLength=int($intergenicLen/2);
		if($halfLength>$INTERGENIC_LENGTH) 
		  {$halfLength=$INTERGENIC_LENGTH}
		my $downStream=substr($intergenicSeq,0,$halfLength);
		my $upstreamBegin=$intergenicLen-$halfLength;
		my $upStream=substr($intergenicSeq,$upstreamBegin,$halfLength);
		my $absoluteUpstream=$intergenicEnd-$halfLength;
		if($upStreamCounter<$N && length($upStream)>=10)
		  {
		    $fastaWriter->addToFasta(">$transcriptId ".
					     "/class=upStreamIntergenic ".
					     "/begin=$absoluteUpstream ".
					     "/len=$halfLength ".
					     "/strand=$strand ",
					     $upStream,
					     $fhUpStreamIntergenic);
		    $upStreamCounter+=length $upStream;
		  }
		if($downStreamCounter<$N && length($downStream)>=10)
		  {
		    $fastaWriter->addToFasta(">$transcriptId ".
					     "/class=downStreamIntergenic ".
					     "/begin=$intergenicBegin ".
					     "/len=$halfLength ".
					     "/strand=$strand ",
					     $downStream,
					     $fhDownStreamIntergenic);
		    $downStreamCounter+=length $downStream;
		  }
	      }
	  }
	undef $transcript;
	undef $codons;
	undef $exons;
      }
  }
close($fhIntrons);
close($fhExons);
close($fhDonors);
close($fhAcceptors);
close($fhUpStreamIntergenic);
close($fhDownStreamIntergenic);
close($fhStartCodons);
close($fhStopCodons);
close($fhDonorFrameShifts);
close($fhAccFrameShifts);
print "done.\n";


#--------------------------------------------------------------
sub attachCodonsToExons
{
  my ($codons)=@_;
  my $numCodons=@$codons;

  # Attach all but the first and last codons to their
  # respective exons
  for(my $i=1 ; $i<$numCodons-1 ; ++$i)
    {
      my $codon=$codons->[$i];
      my $exon=$codon->{exon};
      push @{$exon->{codons}},$codon;
    }
}
#--------------------------------------------------------------
# trimExons(\@exons,$startCodon,$stopCodon,$transcript,$strand)
# Drops exons before start codon or after stop codon.  Also,
# trims portion of exon before start codon, similarly for stop.
sub trimExons
  {
    my ($exons,$startCodon,$stopCodon,$transcript,$strand)=@_;
    my $numExons=@$exons;

    # Adjust transcript begin and end
    if($strand eq "+")
      {
	$transcript->{fivePrime}=$startCodon->{absoluteCoord}+3;
	$transcript->{threePrime}=$stopCodon->{absoluteCoord};
      }
    else # $strand eq "-"
      {
	$transcript->{fivePrime}=$stopCodon->{absoluteCoord};
	$transcript->{threePrime}=$startCodon->{absoluteCoord}-3;
      }

    my $startExon=$startCodon->{exon};
    my $stopExon=$stopCodon->{exon};
    my $strand=$startExon->getStrand();

    # Remove any exons left of the exon containing the start codon
    while($numExons>0)
      {
	last if $exons->[0]==$startExon;
	shift @$exons;
	--$numExons;
      }
    die unless $numExons>0;

    # Trim the first exon so that it begins just after the start codon
    #$startCodon->{relativeCoord}=0;
    my $trimSize=
      ($strand eq "+" ?
       $startCodon->{absoluteCoord}-$startExon->{begin} :
       $startExon->{end}-$startCodon->{absoluteCoord});
    $startExon->trimInitialPortion($trimSize);
    $startExon->trimInitialPortion(3); # remove the start codon

    # Adjust relative coordinates of codons in first exon to account
    # for trimming
    my $codons=$startExon->{codons};
    if($codons) # ATG may have been the only codon in this exon
      {
	my $numCodons=@$codons;
	for(my $i=0 ; $i<$numCodons ; ++$i)
	  {
	    my $codon=$codons->[$i];
	    $codon->{relativeCoord}-=($trimSize+3);
	    if($codon->{relativeCoord}<-3)
	      {die "strand=$strand rel=".$codon->{relativeCoord}.", trimSize=$trimSize, triplet=".$codon->{triplet}}
	  }
      }

    # Remove any exons right of the exon containing the stop codon
    while($numExons>0)
      {
	last if $exons->[$numExons-1]==$stopExon;
	pop @$exons;
	--$numExons;
      }
    die unless $numExons>0;

    # Trim final exon so that it ends just before the stop codon
    my $trimSize=
      ($strand eq "+" ?
       $stopExon->{end}-($stopCodon->{absoluteCoord}+3) :
       $stopCodon->{absoluteCoord}-3-$stopExon->{begin});### why 2?
    $stopExon->trimFinalPortion($trimSize);
    $stopExon->trimFinalPortion(3); # remove the stop codon

    # Renumber the exons to account for deletions
    for(my $i=0 ; $i<$numExons ; ++$i)
      {
	my $exon=$exons->[$i];
	$exon->{order}=$i;
      }
  }
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------



