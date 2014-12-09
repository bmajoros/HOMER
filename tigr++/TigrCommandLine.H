/**************************************************************
TigrCommandLine.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrCommandLine_H
#define INCL_TigrCommandLine_H
#include "TigrVector.H"
#include "TigrString.H"
#include "TigrMap.H"
using namespace std;

class TigrCommandLine
{
public:
  TigrCommandLine(int argc,char *argv[],const char *options);
  TigrString arg(int);
  TigrString optParm(char);
  bool option(char);
  int numArgs();
private:
  int argc;
  char **argv;
  TigrVector<TigrString> args;
  TigrMap<char,bool> usedOption;
  TigrMap<char,TigrString> optionParm;
  bool takesParameter(char c,const char *options);
  void loopThroughOptions(const char *options);
};

#endif
