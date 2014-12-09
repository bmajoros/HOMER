/**************************************************************
TigrPipe.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrPipe_H
#define INCL_TigrPipe_H
#include "TigrFile.H"
using namespace std;

class TigrPipe : public TigrFile
{
public:
  TigrPipe(const TigrString &command,const TigrString &mode);
  virtual void close();
};

#endif
