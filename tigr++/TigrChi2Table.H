/**************************************************************
TigrChi2Table.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrChi2Table_H
#define INCL_TigrChi2Table_H
#include "TigrString.H"
#include "TigrArray1D.H"
#include "TigrArray2D.H"
using namespace std;

static const char *DEFAULT_PATH="";

class TigrChi2Table
{
public:
  TigrChi2Table(TigrString filename=DEFAULT_PATH);
  float lookupChi(int df,float alpha);  
  float lookupP(int df,float chiSquared);  
private:
  void loadFromStream(istream &);

  TigrIntArray1D DFs;
  TigrArray2D<float> table;
  TigrFloatArray1D alphas;
  unsigned alphaToIndex(float alpha);
  unsigned dfToIndex(int df);
};

#endif
