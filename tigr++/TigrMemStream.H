/**************************************************************
TigrMemStream.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrMemStream_H
#define INCL_TigrMemStream_H
#include "TigrString.H"
class TigrMemStream
{
public:
  TigrMemStream(const char *data);
  TigrString readLine();
  bool eof();
private:
  const char *p;
};

#endif
