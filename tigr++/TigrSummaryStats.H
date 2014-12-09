/**************************************************************
TigrSummaryStats.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrSummaryStats_H
#define INCL_TigrSummaryStats_H
#include "TigrVector.H"
using namespace std;

class TigrSummaryStats
{
public:
  TigrSummaryStats(const TigrVector<int> &);
  TigrSummaryStats(const TigrVector<float> &);
  TigrSummaryStats(const TigrVector<double> &);
  double getMax();
  double getMean();
  double getMin();
  double getStdDev();
  int getN();
private:
  int n;
  double mean, stddev, min, max;
  void compute(const TigrVector<double> &);
  void compute(const TigrVector<float> &);
  void compute(const TigrVector<int> &);
};

#endif
