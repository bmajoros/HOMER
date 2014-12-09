/**************************************************************
TigrConfigFile.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrConfigFile_H
#define INCL_TigrConfigFile_H
#include <iostream>
#include "TigrStrTokenizer.H"
#include "TigrMap.H"
using namespace std;

class TigrConfigFile
{
public:
  TigrConfigFile(TigrString fname);
  virtual TigrString lookup(TigrString attr);
  virtual TigrString lookupOrDie(TigrString attr);
  virtual double getDoubleOrDie(TigrString attr);
  virtual float getFloatOrDie(TigrString attr);
  virtual int getIntOrDie(TigrString attr);
  virtual long getLongOrDie(TigrString attr);
  virtual bool getBoolOrDie(TigrString attr);
protected:
  TigrConfigFile();
  void enter(TigrString,TigrString);
private:
  TigrMap<TigrString,TigrString> dict;
  const char *charLookupOrDie(TigrString attr);
  void load(TigrString fname);
  void processLine(TigrStrTokenizer &);
};

#endif
