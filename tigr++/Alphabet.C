/**************************************************************
Alphabet.C
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
using namespace std;
#include "Alphabet.H"
#include <iostream>
#include "TigrString.H"


Alphabet::Alphabet(const char *initializer) : numElements(0)
{
  for(int i=0 ; i<256 ; ++i)
    charToInt[i]=intToChar[i]=static_cast<char>(-1);

  if(initializer)
    {
      int n=strlen(initializer);
      for(int i=0 ; i<n ; ++i) add(initializer[i]);
    }
}



int Alphabet::getNumElements() const
{
  return numElements;
}



int Alphabet::add(char c)
{
  charToInt[static_cast<int>(c)]=static_cast<char>(numElements);
  intToChar[numElements]=c;
  return numElements++;
}



ostream &operator<<(ostream &os,Alphabet &alphabet)
{
  alphabet.printOn(os);
  return os;
}



bool Alphabet::load(istream &is)
{
  char buffer[260];
  is.getline(buffer,260);
  if(is.good())
    {
      char *p=static_cast<char*>(buffer);
      while(*p) add(*p++);
      return true;
    }
  return false;
}



bool Alphabet::save(ostream &os)
{
  for(int i=0 ; i<numElements ; ++i) os << intToChar[i];
  os << endl;
  return os.good() ? true : false;
}



Symbol Alphabet::complement(Symbol s)
{
  return lookup(complement(lookup(s)));
}



char Alphabet::complement(char c)
{
  switch(c)
    {
    case 'A': return 'T';
    case 'G': return 'C';
    case 'C': return 'G';
    case 'T': return 'A';
    case 'N': return 'N';
    }
  throw TigrString("Alphabet::complement:  can't complement base '")+c+"', ASCII "+int(c);
}



void Alphabet::printOn(ostream &os)
{
  for(int i=0 ; i<numElements ; ++i) os << intToChar[i] << ',';
}



bool Alphabet::isDefined(char c) const
{
  return lookup(c)!=-1;
}

