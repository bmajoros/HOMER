/**************************************************************
TigrFastaReader.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrFastaReader_H
#define INCL_TigrFastaReader_H
#include "TigrVector.H"
#include "TigrString.H"
#include "TigrFile.H"
#include "DnaAlphabet.H"
#include "AminoAlphabet.H"
#include "TigrRegex.H"
using namespace std;

/*
  Reads a FASTA file into memory, converts sequence to uppercase,
  and masks strange characters (anything not in the alphabet) using N's.
 */
class TigrFastaReader
{
public:
  TigrFastaReader(const TigrString &filename,Alphabet & =DnaAlphabet::global);
  TigrVector< pair<TigrString,TigrString> > *readAll();
  bool nextSequence(TigrString &defline,TigrString &sequence);
  static void parseDefline(const TigrString &defline,TigrString &id,
			   TigrString &remainder);
  static TigrString getId(const TigrString &defline);
private:
  TigrString cache;
  TigrFile file;
  Alphabet &alphabet;
  static TigrRegex deflineRegex;

  void maskStrangeChars(TigrString &);
};

#endif
