/****************************************************************
 DnaAlphabet.H

 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_DnaAlphabet_H
#define INCL_DnaAlphabet_H
using namespace std;
#include "Alphabet.H"


class DnaAlphabet : public Alphabet
{
public:
  DnaAlphabet();
  static DnaAlphabet global;
private:
};



class LowercaseDnaAlphabet : public Alphabet
{
public:
  LowercaseDnaAlphabet();
  static LowercaseDnaAlphabet global;
private:
};


#endif
