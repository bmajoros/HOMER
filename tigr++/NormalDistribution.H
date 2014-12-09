/****************************************************************
 NormalDistribution.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_NormalDistribution_H
#define INCL_NormalDistribution_H
#include "Constants.H"
#include "ContinuousDistr.H"

class NormalDistribution : public ContinuousDistr
{
  static const double negLogSqrtTwoPi;
  double mu, sigma;
public:
  NormalDistribution(double mu,double sigma);
  virtual double density(double x) const;
  virtual double logDensity(double x) const;
  static inline double density(double x,double mu,double sigma);
  static inline double logDensity(double x,double mu,double sigma);
};



//=======================================================================
//                          inline methods
//=======================================================================
double NormalDistribution::density(double x,double mu,double sigma)
{
  double xMinusMuOverSigma=(x-mu)/sigma;
  return 1/(sqrt(2*PI)*sigma)*exp(-0.5*xMinusMuOverSigma*xMinusMuOverSigma);
}



double NormalDistribution::logDensity(double x,double mu,double sigma)
{
  double xMinusMuOverSigma=(x-mu)/sigma;
  return negLogSqrtTwoPi-log(sigma)-xMinusMuOverSigma*xMinusMuOverSigma/2;
}


#endif
