/**************************************************************
Alphabet.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_Alphabet_H
#define INCL_Alphabet_H
using namespace std;
#include "Symbol.H"
#include <iostream>

class Alphabet
{
public:
  Alphabet(const char * = NULL);
  int add(char);
  bool load(istream &);
  int getNumElements() const;
  inline char lookup(int) const;
  inline int lookup(char) const;
  void printOn(ostream &);
  char complement(char);
  Symbol complement(Symbol);
  bool save(ostream &);
  bool isDefined(char) const;
private:
  char intToChar[256]; // int->char
  char charToInt[256]; // char->int
  int numElements;
};

ostream &operator<<(ostream &,Alphabet &);



char Alphabet::lookup(int index) const
{
  return intToChar[index];
}



int Alphabet::lookup(char c) const
{
  return static_cast<int>(charToInt[static_cast<int>(c)]);
}

#endif
