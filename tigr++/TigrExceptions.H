/**************************************************************
TigrExceptions.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrExceptions_H
#define INCL_TigrExceptions_H
#include "TigrString.H"
using namespace std;

class RootException
{
public:
  RootException(const char *);
  RootException(const TigrString &);
  const TigrString &getMessage() const;
protected:
  TigrString message;
};

class FileErrorException : public RootException
{
public:
  FileErrorException(TigrString filename,TigrString msg)
    : RootException(filename+" "+msg) {}
};

class ArrayIndexException : public RootException
{
public:
  ArrayIndexException(long index,const TigrString &msg);
private:
  void installMessage(long index,const TigrString &);
};

#endif
