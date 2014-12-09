/**************************************************************
TigrGffReader.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrGffReader_H
#define INCL_TigrGffReader_H
#include "TigrVector.H"
#include "TigrString.H"
#include "TigrSet.H"
#include "TigrRegex.H"
#include "TigrGffTranscript.H"
#include "TigrGffFeature.H"
#include "TigrFile.H"
#include "TigrComparator.H"
using namespace std;

class TranscriptComparator : public TigrComparator<TigrGffTranscript*>
{
public:
  bool less(TigrGffTranscript *&,TigrGffTranscript *&);
  bool greater(TigrGffTranscript *&,TigrGffTranscript *&); 
  bool equal(TigrGffTranscript *&,TigrGffTranscript *&);
};


typedef TigrVector<TigrGffTranscript*> TranscriptList;

class TigrGffReader
{
public:
  typedef TigrVector<TigrGffTranscript*> TranscriptList;
  TigrGffReader(const TigrString &filename);
  TigrGffFeature *nextFeature();
  TranscriptList *loadTranscripts();
  TigrMap<TigrString,TranscriptList*> *loadByContig();

private:
  TigrSet<TigrString> exonTypes;
  TigrRegex commentPattern, transgrpRegex;
  TigrFile file;
};

#endif
