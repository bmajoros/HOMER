/**************************************************************
Symbol.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_Symbol_H
#define INCL_Symbol_H
using namespace std;

class Symbol
{
  unsigned char alphabetIndex;
public:
  Symbol(int =0);
  bool operator==(const Symbol &);
  operator int() const;
  int getIndex() const;
  Symbol &operator++();
};


#endif
