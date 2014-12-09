#!/usr/bin/perl
use strict;
use Progress;

my $usage="$0 <*.hmm> <submodels-file>";
die "$usage\n" unless @ARGV==2;
my ($hmmFile,$submodelsFile)=@ARGV;

if(!-e "out") {mkdir "out"}

my $ls=`ls chunks/*.gff`;
my @files=split/\s+/,$ls;
my $numChunks=@files;

my $progress=new Progress;
$progress->start($numChunks);
foreach my $i (1..$numChunks)
  {
    my $chunk="chunks/$i.fasta";
    next unless -e $chunk;
    my $command="homer/homer $hmmFile $submodelsFile $chunk > out/$i.gff";
    print "$command\n";
    system($command);
    my ($timeLeft,$percentDone)=$progress->getProgress($i);

    my $msg="$percentDone\% done.  $timeLeft remaining.";
    my $msgLen=length $msg;
    my $bar='='x$msgLen;
    print "\n$bar\n$msg\n$bar\n\n";
  }
