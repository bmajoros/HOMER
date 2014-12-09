/**************************************************************
TigrTime.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrTime_H
#define INCL_TigrTime_H
#include <sys/time.h>
#include "TigrString.H"
using namespace std;

TigrString getDateAndTime();


class TigrTime
{
public:
  TigrString elapsedTime();
  float elapsedSeconds();
  void startCounting();
  void stopCounting();
private:
  timeval timevalStart, timevalStop;
};

#endif
