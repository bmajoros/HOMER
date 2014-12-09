#!/usr/bin/perl
use strict;

my %blacklist;
my $list=`cd chunks;ls -la`;
my @list=split/\n/,$list;
foreach my $line (@list)
  {
    my @fields=split/\s+/,$line;
    my $size=$fields[4];
    my $name=$fields[8];
    #print "name=[$name] size=[$size]\n";
    next unless $name=~/\.fasta/;
    if($size<300)
      {
	$name=~/(\d+)/;
	my $id=$1;
	$blacklist{$id}=1;
      }
  }

my $browser="/home/bmajoros/genomics/browser/make-figure";
my $base=".";
my $bacsPerPage=4;
my $maxPages=1000;

if(!-d "ps") {system("mkdir ps")}
my $page=1;
for(my $i=1 ; $i<4000 ; $i+=$bacsPerPage, ++$page)
  {
    my $command="$browser 2";
    my $found=0;
    for(my $bac=0 ; $bac<$bacsPerPage ; ++$bac)
      {
	my $id=$i+$bac;
	next if $blacklist{$id};
	my $true="$base/chunks/$id.gff";
	my $ghmm="$base/out/$id.gff";
	next unless -e $ghmm;
	$command.=" $true cDNA\#$id $ghmm Unveil\#$id";
	++$found;
      }
    last unless $found;
    $command.=" ps/page$page.ps";
    system($command);
    last unless $page<$maxPages;
  }

