/**************************************************************
TigrProteinTrans.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrProteinTrans_H
#define INCL_TigrProteinTrans_H
#include "TigrString.H"
#include "TigrMap.H"
#include "Sequence.H"
#include "Alphabet.H"
using namespace std;

class TigrProteinTrans
{
public:
  static TigrString reverseComplement(const TigrString &);
  static TigrString translate(const TigrString &transcript);
  static Sequence *translate(const Sequence &,Alphabet &dnaAlphabet,
			     Alphabet &proteinAlphabet);
  static void translate(const Sequence &source,Alphabet &dnaAlphabet,
			Alphabet &proteinAlphabet,Sequence &dest);
  static char mapCodon(const char *codon);
private:
  static TigrProteinTrans initializer;
  struct mapstrcmp
  {
    bool operator()(const char *p,const char *q)
      {	return strcmp(p,q)<0; }
  };
  static TigrMap<const char *,char,mapstrcmp> codonMap;
  TigrProteinTrans();
  static char complement(char);
  void setupCodonMap();
};

#endif
