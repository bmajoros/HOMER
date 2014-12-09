/**************************************************************
TigrGffTranscript.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrGffTranscript_H
#define INCL_TigrGffTranscript_H
#include <iostream>
#include "TigrVector.H"
#include "TigrString.H"
#include "TigrGffExon.H"
using namespace std;

class TigrGffTranscript
{
public:
  TigrGffTranscript(const TigrString &transcriptId,
		    const TigrString &substrate,
		char strand,const TigrString &source);
  virtual ~TigrGffTranscript();
  TigrGffExon &getIthExon(int);
  TigrGffFeature *getStartCodon();
  TigrGffFeature *getStopCodon();
  char getStrand();
  const TigrString &getSource() const;
  const TigrString &getSubstrate() const;
  const TigrString &getTranscriptId() const;
  double getScore() const;
  int getBegin() const;
  int getEnd() const;
  int getNumExons() const;
  void addExon(TigrGffExon *);
  void printOn(ostream &) const;
  void setScore(double);
  void setStartCodon(TigrGffFeature *);
  void setStopCodon(TigrGffFeature *);
  void setStrand(char);
  void sortExons();
  void setExonTypes();
  void toGff(ostream &);
  TigrVector<TigrGffExon*>::iterator getExons() {return exons.begin();}
  TigrVector<TigrGffExon*>::iterator getExonsEnd() {return exons.end();}
private:
  int begin, end;  
  double score;
  char strand;
  bool hasScore;
  TigrVector<TigrGffExon*> exons;
  TigrString transcriptId;
  TigrString substrate;    
  TigrString source;       
  TigrGffFeature *startCodon, *stopCodon;
};

ostream &operator<<(ostream &,const TigrGffTranscript &);

#endif
