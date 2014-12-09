/****************************************************************
 AminoAlphabet.C

 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
using namespace std;
#include "AminoAlphabet.H"
#include <iostream>


AminoAlphabet AminoAlphabet::global;


AminoAlphabet::AminoAlphabet()
  : Alphabet("*ARNDCQEGHILKMFPSTWYV") // plus BZX ...
{
}

