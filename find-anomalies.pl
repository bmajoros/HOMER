#!/usr/bin/perl
use strict;

my $biggest;
for(my $seed=1 ; $seed<500 ; ++$seed)
  {
    system("/home/bmajoros/GF/homer/baum-welch intron.hmms introns.fasta 25 junk 1 3 -s $seed >& tmp.1");
    open(PIPE,"tmp.1") || die;
    while(<PIPE>)
      {
	if(/>>> (\S+) <<</)
	  {
	    my $decrease=abs($1);
	    if(!defined($biggest) || $decrease>$biggest) 
	      {$biggest=$decrease}
	    print "$decrease   (biggest so far: $biggest)  seed=$seed\n";
	  }
      }
    close(PIPE);
  }




