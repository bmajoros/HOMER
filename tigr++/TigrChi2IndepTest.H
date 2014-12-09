/**************************************************************
TigrChi2IndepTest.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrChi2IndepTest_H
#define INCL_TigrChi2IndepTest_H
using namespace std;
#include "TigrArray2D.H"
#include "TigrChi2Table.H"

class TigrChi2IndepTest
{
public:
  TigrChi2IndepTest(TigrArray2D<int> &contingencyTable,TigrChi2Table &);
  bool areIndependent();
  double getChiSquared();
  double getP();
  inline bool areDependent() {return !areIndependent();}
private:
  double P, chiSquared;
  bool independent;
  void runTest(TigrArray2D<int> &,TigrChi2Table &);
};

#endif
