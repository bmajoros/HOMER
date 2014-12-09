/**************************************************************
TigrMultRegress.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrMultRegress_H
#define INCL_TigrMultRegress_H
#include "TigrDblMatrix.H"
using namespace std;

class TigrMultRegress
{
public:
  bool regress(const TigrDblMatrix &Y,const TigrDblMatrix &X,
	       TigrDblMatrix &coef);
private:
  bool inversionRegression(const TigrDblMatrix &Y,const TigrDblMatrix &X,
			   TigrDblMatrix &coef);
};

#endif
