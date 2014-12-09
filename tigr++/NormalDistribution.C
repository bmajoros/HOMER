/****************************************************************
 NormalDistribution.C
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#include "NormalDistribution.H"
#include <iostream>


const double NormalDistribution::negLogSqrtTwoPi=-log((long double)2*PI)/2;


NormalDistribution::NormalDistribution(double mu,double sigma)
  : mu(mu), sigma(sigma)
{
}



double NormalDistribution::density(double x) const
{
  return density(x,mu,sigma);
}



double NormalDistribution::logDensity(double x) const
{
  return logDensity(x,mu,sigma);
}


