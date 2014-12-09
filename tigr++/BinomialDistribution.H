/****************************************************************
 BinomialDistribution.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_BinomialDistribution_H
#define INCL_BinomialDistribution_H
#include <math.h>

class BinomialDistribution
{
public:
  static inline double density(int exactNumSuccesses,int numTrials,
			       double probSuccess);
  static double rightTailedPValue(int thisManyOrMoreSuccesses,int numTrials,
				  double probabilityOfSuccess);
};



//=======================================================================
//                          inline methods
//=======================================================================
double BinomialDistribution::density(int x,int n,double p)
{
  // term 1: sum(d=x+1 to n) log d
  double term1=0;
  for(int d=x+1 ; d<=n ; ++d) 
    term1+=log(float(d));
  
  // term 2: -sum(d=1 to n-x) log d
  double term2=0;
  for(int d=1 ; d<=n-x ; ++d) 
    term1-=log(float(d));
  
  // term 3: x log p
  double term3=x*log(p);
  
  // term 4: (n-x) log (1-p)
  double term4=(n-x)*log(1-p);
  
  // probability = e to the (term1+term2+term3+term4)
  return exp(term1+term2+term3+term4);
}


#endif
