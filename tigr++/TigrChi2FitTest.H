/**************************************************************
TigrChi2FitTest.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrChi2FitTest_H
#define INCL_TigrChi2FitTest_H
using namespace std;
#include "TigrVector.H"
#include "TigrChi2Table.H"

class TigrChi2FitTest
{
public:
  TigrChi2FitTest(TigrVector<int> &observedCounts,
		    TigrVector<int> &expectedCounts,
		    TigrChi2Table &);
  bool goodFit();
  double getChiSquared();
  double getP();
private:
  double P, chiSquared;
  bool fit;
  void performTest(TigrVector<int> &observedCounts,
		   TigrVector<int> &expectedCounts,
		   TigrChi2Table &);
};

#endif
