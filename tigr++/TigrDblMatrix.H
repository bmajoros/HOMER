/**************************************************************
TigrDblMatrix.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrDblMatrix_H
#define INCL_TigrDblMatrix_H
#include <iostream>
#include "TigrArray2D.H"
#include "TigrArray1D.H"
using namespace std;
class TigrDblMatrix;

void multiply(TigrDblMatrix &leftM,TigrDblMatrix &rightM,
	      TigrDblMatrix &resultMatrix);

class TigrDblMatrix
{
public:
  TigrDblMatrix(int rows,int columns);
  TigrDblMatrix(const TigrDblMatrix &other);
  virtual double operator()(int row,int col) const;
  virtual double &operator()(int row,int col);
  TigrDblMatrix &operator=(TigrDblMatrix &);
  bool invert(TigrDblMatrix &resultingMatrix) const;
  virtual int getNumColumns() const;
  virtual int getNumRows() const;
  void addRowMultiple(int sourceRow,int destinationRow,double factor);
  void getColumn(int column,TigrDblArray1D &into) const;
  void getRow(int row,TigrDblArray1D &into) const;
  void multiplyRowBy(int whichRow,double factor);
  void printOn(ostream &);
  void setAllTo(double d);
  void swapRows(int r1,int r2);
  void times(const TigrDblMatrix &,TigrDblMatrix &resultMatrix) const;
  void transpose(TigrDblMatrix &resultMatrix) const;
private:
  TigrArray2D<double> theArray;
};

ostream &operator<<(ostream &,TigrDblMatrix &);

#endif
