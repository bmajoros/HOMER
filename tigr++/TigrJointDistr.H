/**************************************************************
TigrJointDistr.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrJointDistr_H
#define INCL_TigrJointDistr_H
#include "TigrArray2D.H"
using namespace std;

class TigrJointDistr : public TigrArray2D<double>
{
public:
  TigrJointDistr(int X,int Y);
  double getMarginalX(int x);
  double getMarginalY(int y);
  void computeMarginals();
private:
  TigrDblArray1D marginalX, marginalY;
};

#endif
