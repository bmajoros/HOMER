/**************************************************************
TigrMutInfo.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrMutInfo_H
#define INCL_TigrMutInfo_H
#include "TigrJointDistr.H"
#include "TigrContingencyTbl.H"
using namespace std;

class TigrMutInfo
{
public:
  static double compute(const TigrContingencyTbl &);
  static double compute(const TigrJointDistr &);
};

#endif
