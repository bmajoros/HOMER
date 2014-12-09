#!/usr/bin/perl
use strict;
my $unveil="unveil";
use Getopt::Std;

our($opt_i);
getopts("i:");

my $EM_ITERATIONS=5;
if($opt_i) {$EM_ITERATIONS=$opt_i}

determ("exon.hmms","exons.fasta","exon.hmm");
determ("splice-junction.hmms","donors.fasta","donor.hmm");
determ("splice-junction.hmms","acceptors.fasta","acceptor.hmm");
determ("start-codon.hmms","start-codons.fasta","start-codon.hmm");
baumWelch("frame-shift.hmms","donor-frame-shifts.fasta",$EM_ITERATIONS,
	  "donor-frame-shift.hmm");
baumWelch("frame-shift.hmms","acceptor-frame-shifts.fasta",$EM_ITERATIONS,
	  "acceptor-frame-shift.hmm");
baumWelch("intron.hmms","introns.fasta",$EM_ITERATIONS,"intron.hmm");
baumWelch("stop-codon.hmms","stop-codons.fasta",$EM_ITERATIONS,
	  "stop-codon.hmm");
baumWelch("upstream-intergenic.hmms","upstream-intergenic.fasta",
	  $EM_ITERATIONS,"upstream-intergenic.hmm");
baumWelch("downstream-intergenic.hmms","downstream-intergenic.fasta",
	  $EM_ITERATIONS,"downstream-intergenic.hmm");


#---------------------------------------------------------------------
sub baumWelch
  {
    my ($structFile,$exampleFile,$numIterations,$outfile)=@_;
    system("date");
    print STDERR "Running baum-welch on $exampleFile...\n";
    system("$unveil/baum-welch $structFile $exampleFile $numIterations $outfile");
  }
#---------------------------------------------------------------------
sub determ
  {
    my ($structFile,$exampleFile,$outfile)=@_;
    system("date");
    print STDERR "Running deterministic-trainer on $exampleFile...\n";
    system("$unveil/deterministic-trainer $structFile $exampleFile $outfile");
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
