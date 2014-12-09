/**************************************************************
TigrGffExon.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrGffExon_H
#define INCL_TigrGffExon_H
#include "TigrMap.H"
#include "TigrGffFeature.H"
using namespace std;

enum ExonType
  {
    ET_EXON,
    ET_INITIAL_EXON,
    ET_INTERNAL_EXON,
    ET_FINAL_EXON,
    ET_SINGLE_EXON
  };
class TigrGffTranscript;

ostream &operator<<(ostream &,ExonType);

class ExonTypeInitializer
{
private:
  static ExonTypeInitializer e;
  ExonTypeInitializer();
};

class TigrGffExon
{
public:
  TigrGffExon(TigrGffFeature &,TigrGffTranscript &);
  TigrGffExon(ExonType,int begin,int end,TigrGffTranscript &parent,
	  bool hasScore,double score,bool hasFrame,int frame);
  ExonType getExonType() const;
  void changeExonType(ExonType);
  TigrGffTranscript &getParent();
  char getStrand() const;
  const TigrString &getSource() const;
  const TigrString &getSubstrate() const;
  int getBegin();  
  int getEnd();    
  void toGff(ostream &);
private:
  static TigrMap<TigrString,ExonType> exonTypeNames;
  int frame;  
  int begin, end;            
  double score;
  bool hasScore;  
  bool hasFrame;  
  TigrGffTranscript &parent;
  ExonType exonType;
  friend class ExonTypeInitializer;
  static void initExonTypeNames();
};

#endif
