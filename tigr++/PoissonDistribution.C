/****************************************************************
 PoissonDistribution.C
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#include "PoissonDistribution.H"
#include <iostream>


PoissonDistribution::PoissonDistribution(int lambda)
  : lambda(lambda)
{
}



double PoissonDistribution::density(int x)
{
  return density(x,lambda);
}



double PoissonDistribution::logDensity(int x)
{
  return logDensity(x,lambda);
}


