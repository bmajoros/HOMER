/**************************************************************
TigrLinearFunc.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrLinearFunc_H
#define INCL_TigrLinearFunc_H
using namespace std;

class TigrLinearFunc
{
public:
  TigrLinearFunc(double slope,double intercept,
		 double coefOfDeterm);
  TigrLinearFunc();
  double operator()(double) const;
  TigrLinearFunc &operator=(const TigrLinearFunc &);
  double getCoefDeterm() const;
private:
  double slope, intercept;
  double coefOfDeterm;
};

#endif
