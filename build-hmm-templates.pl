#!/usr/bin/perl
use strict;

##############################################################
# submodels.txt
##############################################################
my $dir=`pwd`;
chop $dir;
open(OUT,">submodels.txt") || die;
print OUT "1[upstream-intergenic]=$dir/upstream-intergenic.hmm
2[downstream-intergenic]=$dir/downstream-intergenic.hmm
3[start]=$dir/start-codon.hmm
4[stop]=$dir/stop-codon.hmm
5[exon]=$dir/exon.hmm
6[frameshift]=$dir/donor-frame-shift.hmm
7[frameshift]=$dir/acceptor-frame-shift.hmm
8[donor-2]=$dir/donor.hmm
9[acceptor-1]=$dir/acceptor.hmm
10=$dir/intron.hmm
";
close(OUT);

##############################################################
# metamodel.hmms
##############################################################
open(OUT,">metamodel.hmms") || die;
print OUT "0 -> 1 : 1
1 -> 3 : 1
2 -> 1 : 0.9
2 -> 0 : 0.1
3 -> 5 : 1
4 -> 2 : 1
5 -> 6 : 0.6
5 -> 4 : 0.4
6 -> 8 : 1
7 -> 5 : 1
8 -> 10 : 1
9 -> 7 : 1
10 -> 9 : 1
";
close(OUT);

##############################################################
# Start codon
##############################################################
open(OUT,">start-codon.hmms");
print OUT "
0 -> 1 : 1.0
1 -> 2 : 1.0
2 -> 3 : 1.0
3 -> 0 : 1.0
state 1 : A=0.97 T=0.01 C=0.01 G=0.01
state 2 : T=0.97 A=0.01 C=0.01 G=0.01
state 3 : G=0.97 A=0.01 T=0.01 C=0.01
";
close(OUT);

##############################################################
# Stop codon
##############################################################
open(OUT,">stop-codon.hmms");
print OUT "
0 -> 1 : 1
1 -> 2 : 0.5
1 -> 5 : 0.5
2 -> 3 : 0.5
2 -> 4 : 0.5
3 -> 0 : 1
4 -> 0 : 1
5 -> 6 : 1
6 -> 0 : 1
state 1 : T=0.097 A=0.001 G=0.001 C=0.001
state 2 : A=0.097 T=0.001 G=0.001 C=0.001
state 3 : A=0.097 T=0.001 G=0.001 C=0.001
state 4 : G=0.097 A=0.001 T=0.001 C=0.001
state 5 : G=0.097 A=0.001 T=0.001 C=0.001
state 6 : A=0.097 T=0.001 G=0.001 C=0.001
";
close(OUT);

##############################################################
# Generic frame shift
##############################################################
open(OUT,">frame-shift.hmms");
print OUT "
0 -> 0 : 0.5
0 -> 1 : 0.5
1 -> 0 : 0.5
1 -> 2 : 0.5
2 -> 0 : 1
";
close(OUT);

##############################################################
# Generic splice junction
##############################################################
open(OUT,">splice-junction.hmms");
print OUT "
0 -> 1 : 1
1 -> 2 : 1
2 -> 3 : 1
3 -> 4 : 1
4 -> 5 : 1
5 -> 6 : 1
6 -> 7 : 1
7 -> 0 : 1
";
close(OUT);

##############################################################
# Intron model
##############################################################
open(OUT,">intron.hmms");
print OUT "
0 -> 1 : 0.33333
0 -> 2 : 0.33333
0 -> 3 : 0.33334
1 -> 4 : 1
2 -> 1 : 1
3 -> 2 : 1
4 -> 5 : 1
5 -> 0 : 0.25
5 -> 6 : 0.25
5 -> 9 : 0.25
5 -> 10 : 0.25
6 -> 7 : 1
7 -> 8 : 1
8 -> 0 : 0.33333
8 -> 1 : 0.33333
8 -> 9 : 0.33334
9 -> 10 : 1
10 -> 0 : 1
";
close(OUT);

##############################################################
# Exon
##############################################################
open(OUT,">exon.hmms");
my @alphabet=('A','T','C','G');
my (%trans,%emit);
my $nextStateId=1;
for(my $i=0 ; $i<4 ; ++$i)
  {
    my $firstSymbol=$alphabet[$i];
    my $firstState=$nextStateId++;
    $emit{$firstState}=$firstSymbol;
    push @{$trans{0}},$firstState;
    for(my $j=0 ; $j<4 ; ++$j)
      {
	my $secondSymbol=$alphabet[$j];
	my $secondState=$nextStateId++;
	$emit{$secondState}=$secondSymbol;
	push @{$trans{$firstState}},$secondState;

	for(my $k=0 ; $k<4 ; ++$k)
	  {
	    my $thirdSymbol=$alphabet[$k];
	    my $thirdState=$nextStateId++;
	    $emit{$thirdState}=$thirdSymbol;
	    push @{$trans{$secondState}},$thirdState;
	    push @{$trans{$thirdState}},0;
	    push @{$trans{$thirdState}},1;
	    push @{$trans{$thirdState}},22;
	    push @{$trans{$thirdState}},43;
	    push @{$trans{$thirdState}},64;
	  }
      }
  }
for(my $i=0 ; $i<$nextStateId ; ++$i)
  {
    my $transitions=$trans{$i};
    my $numTo=@$transitions;
    my $p=1/$numTo;
    foreach my $to (@$transitions)
      {
	print OUT "$i -> $to : $p\n";
      }
  }
for(my $i=1 ; $i<$nextStateId ; ++$i)
  {
    my $symbol=$emit{$i};
    print OUT "state $i : $symbol=1 ";
    foreach my $other (@alphabet)
      {
	next if $other eq $symbol;
	print OUT "$other=0 ";
      }
    print OUT "\n";
  }
close(OUT);

##############################################################
# Downstream intergenic model
##############################################################
open(OUT,">downstream-intergenic.hmms");
print OUT "
0 -> 1 : 1
1 -> 2 : 1
2 -> 3 : 1
3 -> 4 : 1
4 -> 5 : 1
5 -> 6 : 1
6 -> 7 : 1
7 -> 8 : 1
8 -> 9 : 1
9 -> 10 : 1
10 -> 10 : 0.5
10 -> 0 : 0.5
";
close(OUT);

##############################################################
# Upstream intergenic model
##############################################################
open(OUT,">upstream-intergenic.hmms");
print OUT "
0 -> 1 : 1
1 -> 1 : 0.5
1 -> 2 : 0.5
2 -> 3 : 1
3 -> 4 : 1
4 -> 5 : 1
5 -> 6 : 1
6 -> 7 : 1
7 -> 8 : 1
8 -> 9 : 1
9 -> 10 : 1
10 -> 0 : 1
";
close(OUT);
