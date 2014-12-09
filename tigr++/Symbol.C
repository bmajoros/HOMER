/**************************************************************
Symbol.C
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/

using namespace std;
#include "Symbol.H"
#include <iostream>
#include "TigrString.H"


Symbol::Symbol(int alphabetIndex) : alphabetIndex(alphabetIndex)
{
  if(alphabetIndex>255 || alphabetIndex<0) 
    throw TigrString("")+alphabetIndex+
      " is too large for char in Symbol::Symbol";
}



bool Symbol::operator==(const Symbol &other) 
{ 
  return alphabetIndex==other.alphabetIndex;
}



Symbol::operator int() const 
{ 
  return alphabetIndex;
}



int Symbol::getIndex() const
{
  return alphabetIndex;
}



Symbol &Symbol::operator++() 
{ 
  ++alphabetIndex; 
  return *this;
}



