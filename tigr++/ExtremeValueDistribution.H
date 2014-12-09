/****************************************************************
 ExtremeValueDistribution.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ExtremeValueDistribution_H
#define INCL_ExtremeValueDistribution_H
#include <iostream>
#include <math.h>
#include "ContinuousDistr.H"
using namespace std;

class ExtremeValueDistribution : public ContinuousDistr
{
  double alpha, beta;
public:
  ExtremeValueDistribution(double alpha,double beta);
  virtual double density(double x) const;
  virtual double logDensity(double x) const;
  static inline double density(double x,double alpha,double beta);
  static inline double logDensity(double x,double alpha,double beta);
};



//=======================================================================
//                          inline methods
//=======================================================================
double ExtremeValueDistribution::density(double x,double alpha,double beta)
{
  double aMinusXOverB=(alpha-x)/beta;
  return exp(aMinusXOverB-exp(aMinusXOverB))/beta;
}



double ExtremeValueDistribution::logDensity(double x,double alpha,double beta)
{
  double aMinusXOverB=(alpha-x)/beta;
  return aMinusXOverB-exp(aMinusXOverB)-log(beta);
}


#endif
