/****************************************************************
 ExtremeValueDistribution.C
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#include "ExtremeValueDistribution.H"
#include <iostream>


ExtremeValueDistribution::ExtremeValueDistribution(double alpha,double beta)
  : alpha(alpha), beta(beta)
{
}



double ExtremeValueDistribution::density(double x) const
{
  return density(x,alpha,beta);
}



double ExtremeValueDistribution::logDensity(double x) const
{
  return logDensity(x,alpha,beta);
}



