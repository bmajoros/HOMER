#!/usr/bin/perl
use strict;

my $ls=`ls chunks`;
my @files=split/\s+/,$ls;
my ($codingLength,$totalLength);
foreach my $file (@files)
  {
    next unless $file=~/(\S+)\.gff$/;
    my $filestem=$1;
    my $gff=$file;
    my $fasta="chunks/$filestem.fasta";
    my $fastaLen=`fasta-seq-length.pl $fasta`;
    $fastaLen=~/(\d+)/;
    $fastaLen=$1;
    my $stats=`gff-stats.pl chunks/$gff`;
    $stats=~/transcript length:\s*(\d+)/ 
      || die "can't parse stats: $stats\n";
    my $transcriptLen=$1;
    $codingLength+=$transcriptLen;
    $totalLength+=$fastaLen;
  }
my $p=$codingLength/$totalLength;
my $density=int(100*$p+0.5);
print "coding density=$density\%\n";
my $baseline=int(100*(2*$p*$p-2*$p+1)+0.5);
print "baseline=$baseline\%\n";



