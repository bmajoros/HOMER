/**************************************************************
TigrEntropy.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrEntropy_H
#define INCL_TigrEntropy_H
#include "TigrVector.H"
using namespace std;

class TigrEntropy
{
public:
  static double entropy(TigrVector<double> &);
  static double relativeEntropy(TigrVector<double> &,TigrVector<double> &);
  static double crossEntropy(TigrVector<double> &,TigrVector<double> &);
};

#endif
