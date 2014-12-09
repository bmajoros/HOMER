/**************************************************************
TigrStacktrace.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrStacktrace_H
#define INCL_TigrStacktrace_H
#include "TigrString.H"
using namespace std;

class Crasher  
{
public:
  virtual void crash() {};
};

class TigrStacktrace  
{
public:
  TigrStacktrace(const char *message);
  TigrStacktrace(const TigrString &);
};

#endif
