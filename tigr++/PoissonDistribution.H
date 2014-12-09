/****************************************************************
 PoissonDistribution.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_PoissonDistribution_H
#define INCL_PoissonDistribution_H
#include "Factorial.H"


class PoissonDistribution
{
  int lambda;
public:
  PoissonDistribution(int expected);
  double density(int observed);
  double logDensity(int observed);
  static inline double density(int observed,int expected);
  static inline double logDensity(int observed,int expected);
};



//=======================================================================
//                          inline methods
//=======================================================================
double PoissonDistribution::density(int x,int lambda)
{
  return exp(logDensity(x,lambda));
}



double PoissonDistribution::logDensity(int x,int lambda)
{
  return -lambda+x*log((long double)lambda)-Factorial::global.logFactorial(x);
}


#endif
