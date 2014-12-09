/**************************************************************
TigrLinRegressor.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrLinRegressor_H
#define INCL_TigrLinRegressor_H
#include <utility>
#include "TigrVector.H"
#include "TigrLinearFunc.H"
using namespace std;

typedef pair<double,double> DblPoint;
typedef TigrVector<DblPoint> DblPointVector;

class TigrLinRegressor
{
public:
  TigrLinearFunc regress(DblPointVector &);
  TigrLinearFunc regress(TigrVector<double> &);
private:
  TigrLinearFunc compute(double sumX,double sumXX,double sumY,
			 double sumYY,double sumXY,int n);
};

#endif
