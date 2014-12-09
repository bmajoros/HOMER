/****************************************************************
 BinomialDistribution.C
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#include "BinomialDistribution.H"
#include <iostream>


double BinomialDistribution::rightTailedPValue(int successes,int n,
					       double p)
{
  /*
    Uses binomial distribution to compute probability of seeing
    #successes in n trials (or more successes in n trials; hence
    the summation), assuming trials are independent and identically
    distributed (iid).  Thus, this function provides the P-value for
    a right-tailed test.
   */

  double P=0;
  for(int x=successes ; x<=n ; ++x)
    P+=density(x,n,p);

  return P;
}


