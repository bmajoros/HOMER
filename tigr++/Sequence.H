/**************************************************************
Sequence.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_Sequence_H
#define INCL_Sequence_H
using namespace std;
#include <iostream>
#include <ctype.h>
#include "TigrString.H"
#include "TigrVector.H"
#include "Alphabet.H"
#include "Symbol.H"
#include "TigrRegex.H"


class Sequence
{
public:
  static Sequence *load(const TigrString &filename,Alphabet &,
			TigrString &substrateId);
  static Sequence *load(const TigrString &filename,Alphabet &);
  Sequence(const TigrString &,Alphabet &);
  Sequence();
  Sequence(const char *,Alphabet &);
  virtual ~Sequence() {}
  void copyFrom(const TigrString &,Alphabet &);
  void setPhase(int);
  int getPhase() const;
  void append(Symbol);
  void append(const Sequence &);
  void prepend(Symbol);
  Symbol operator [](int) const;
  Symbol &operator [](int);
  const Sequence &operator=(const Sequence &);
  int getLength() const;
  void printOn(ostream &) const;
  void printOn(ostream &,Alphabet &) const;
  TigrString *toString(Alphabet &,int startingAt=0) const;
  void clear();
  virtual void getSubsequence(int begin,int len,Sequence &) const;
  bool subsequenceOccursAt(const Sequence &,int at) const;
  Sequence *reverseComplement(Alphabet &) const;
  void reverseComplement(Alphabet &,Sequence &) const;
  bool isEmpty() const {return symbols.isEmpty();}
  int countOccurrences(Symbol);
private:
  static TigrRegex deflineRegex;
  void append(Alphabet &,const char *);
  TigrVector<Symbol> symbols;
  int phase;
};

ostream &operator<<(ostream &,const Sequence &);

#endif
