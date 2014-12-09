#!/usr/bin/perl
use strict;
use GffTranscriptReader;

die unless -e "chunks";
my $gffReader=new GffTranscriptReader;
my $ls=`ls chunks`;
my @files=split/\s+/,$ls;
foreach my $file (@files)
  {
    next unless($file=~/(\S+)\.gff/);
    my $filestem=$1;
    my $gffFile="chunks/$filestem.gff";
    my $fastaFile="chunks/$filestem.fasta";
    my $transcripts=$gffReader->loadGFF($gffFile);
    my $transcript=$transcripts->[0];
    my $strand=$transcript->getStrand();
    next unless $strand eq "-";
    my $seqLen=0+`fasta-seq-length.pl $fastaFile`;
    system("revcomp-gff.pl $gffFile $seqLen > tmp.1 ; mv tmp.1 $gffFile");
    system("revcomp-fasta.pl $fastaFile > tmp.1 ; mv tmp.1 $fastaFile");
  }


