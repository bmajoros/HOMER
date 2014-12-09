/****************************************************************
 TigrSequenceEntropy.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_TigrSequenceEntropy_H
#define INCL_TigrSequenceEntropy_H

#include "Sequence.H"
#include "TigrString.H"

class TigrSequenceEntropy
{
public:
  static double entropy(const Sequence &,int order,double &maxEntropy);
  static double entropy(const TigrString &,int order,double &maxEntropy);
};


#endif
