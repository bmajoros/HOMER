/****************************************************************
 ContinuousDistr.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ContinuousDistr_H
#define INCL_ContinuousDistr_H

class ContinuousDistr
{
public:
  virtual double density(double x) const=0;
  virtual double logDensity(double x) const=0;
};


#endif
