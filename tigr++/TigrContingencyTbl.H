/**************************************************************
TigrContingencyTbl.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrContingencyTbl_H
#define INCL_TigrContingencyTbl_H
#include "TigrArray2D.H"
#include "TigrArray1D.H"
using namespace std;

class TigrContingencyTbl : public TigrIntArray2D
{
public:
  TigrContingencyTbl(int width,int height);
  int getColumnTotal(int column);
  int getGrandTotal();
  int getRowTotal(int row);
  void computeTotals();
private:
  int grandTotal;
  TigrArray1D<int> columnTotals, rowTotals;
};

#endif
