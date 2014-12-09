/**************************************************************
TigrRealVector.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrRealVector_H
#define INCL_TigrRealVector_H
#include "TigrDblMatrix.H"
using namespace std;

class TigrRealVector
{
public:
  TigrRealVector(int dimension);
  double operator[](int index) const;
  double &operator[](int index);
  double distanceTo(const TigrRealVector &other);
  double dotProduct(const TigrRealVector &with);
  double getAngle(const TigrRealVector &from);
  double getMagnitude();
  inline int getDimension();
  void scale(double factor);
  void unitize();
private:
  TigrDblMatrix theMatrix;
};

#endif
