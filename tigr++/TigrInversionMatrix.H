/**************************************************************
TigrInversionMatrix.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrInversionMatrix_H
#define INCL_TigrInversionMatrix_H
#include "TigrDblMatrix.H"
using namespace std;

class TigrInversionMatrix : public TigrDblMatrix
{
public:
  TigrInversionMatrix(int rows);
  virtual double &operator()(int row,int col);
  bool detectNonInvertible();
  virtual int getNumColumns();
  void getInverted(TigrDblMatrix &);
  void install(TigrDblMatrix &);
private:
  void setupIdentMatrix();
};

#endif
