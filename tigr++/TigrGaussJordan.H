/**************************************************************
TigrGaussJordan.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrGaussJordan_H
#define INCL_TigrGaussJordan_H
#include "TigrInversionMatrix.H"
using namespace std;

class TigrGaussJordan
{
public:
  static bool invert(TigrDblMatrix &thisMatrix,TigrDblMatrix &intoHere);
private:
  static int getNonzeroRow(int,TigrInversionMatrix &);
  static void zeroOut(int,TigrInversionMatrix &);
};

#endif
