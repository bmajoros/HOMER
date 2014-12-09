/**************************************************************
TigrGffFeature.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrGffFeature_H
#define INCL_TigrGffFeature_H
#include "TigrVector.H"
#include "TigrString.H"
using namespace std;

class TigrGffFeature
{
public:
  TigrGffFeature(const TigrString &rawLine,const TigrString &substrate,
	     const TigrString &source,const TigrString &featureType,
	     int begin,int end,double score,bool hasScore,char strand,
	     int frame,bool hasFrame);
  TigrGffFeature(const TigrString &rawLine);
  TigrString toGff() const;
  TigrVector<TigrString> &getExtraFields();
  bool hasExtraFields() const;
  bool isFramed() const;  
  bool isScored() const;  
  bool isStranded() const;  
  char getStrand() const;  
  const TigrString &getFeatureType() const;
  const TigrString &getRawField(int index);
  const TigrString &getRawLine() const;
  const TigrString &getSource() const;
  const TigrString &getSubstrate() const;
  double getScore() const;
  int getBegin() const;  // 0-based, always less than end
  int getEnd() const;    // 0-based, always greater than begin
  int getFrame() const;  
private:
  int frame;  
  int begin, end;  
  double score;
  char strand;  
  bool hasScore;  
  bool hasFrame;  
  TigrVector<TigrString> extraFields, allFields;
  TigrString substrate;  
  TigrString source;  
  TigrString rawLine;  
  TigrString featureType;  
  void parseLine(const TigrString &);
};

ostream &operator<<(ostream&,TigrGffFeature&);

#endif
