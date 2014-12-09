/****************************************************************
 AminoAlphabet.H

 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_AminoAlphabet_H
#define INCL_AminoAlphabet_H
using namespace std;
#include "Alphabet.H"

class AminoAlphabet : public Alphabet
{
public:
  AminoAlphabet();
  static AminoAlphabet global;
private:
};


#endif
