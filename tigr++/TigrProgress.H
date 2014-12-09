/**************************************************************
TigrProgress.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrProgress_H
#define INCL_TigrProgress_H
#include "TigrString.H"
using namespace std;

class TigrProgress
{
public:
  TigrProgress();
  TigrString getProgress(unsigned long workDone);
  void start(unsigned long totalWork);
private:
  unsigned long totalWork;
  long startTime;
};

#endif
