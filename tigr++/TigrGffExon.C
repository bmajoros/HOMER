#include <iostream>
#include "TigrGffTranscript.H"
#include "TigrGffExon.H"
using namespace std;



TigrMap<TigrString,ExonType> TigrGffExon::exonTypeNames;
ExonTypeInitializer ExonTypeInitializer::e;


TigrGffExon::TigrGffExon(TigrGffFeature &feature,TigrGffTranscript &parent)
  : parent(parent)
{
  exonType=exonTypeNames[feature.getFeatureType()];
  begin=feature.getBegin();
  end=feature.getEnd();
  score=feature.getScore();
  hasScore=feature.isScored();
  frame=feature.getFrame();
  hasFrame=feature.isFramed();
}



TigrGffExon::TigrGffExon(ExonType type,int begin,int end,
			 TigrGffTranscript &parent,bool hasScore,
			 double score,bool hasFrame,int frame)
  : exonType(type), begin(begin), end(end), parent(parent), 
    hasScore(hasScore), score(score), hasFrame(hasFrame), frame(frame)
{
}



ostream &operator<<(ostream &os,ExonType t)
{
  switch(t)
    {
    case ET_EXON:          os<<"exon";            break;
    case ET_INITIAL_EXON:  os<<"initial-exon";    break;
    case ET_INTERNAL_EXON: os<<"internal-exon";   break;
    case ET_FINAL_EXON:    os<<"final-exon";      break;
    case ET_SINGLE_EXON:   os<<"single-exon";     break;
    }
  return os;
}



ExonType TigrGffExon::getExonType() const
{
  return exonType;
}



ExonTypeInitializer::ExonTypeInitializer() 
{
  TigrGffExon::initExonTypeNames();
}



TigrGffTranscript &TigrGffExon::getParent()
{
  return parent;
}



char TigrGffExon::getStrand() const
{
  return parent.getStrand();
}



const TigrString &TigrGffExon::getSource() const
{
  return parent.getSource();
}



const TigrString &TigrGffExon::getSubstrate() const
{
  return parent.getSubstrate();
}



int TigrGffExon::getBegin()  
{
  return begin;
}



int TigrGffExon::getEnd()  
{
  return end;
}



void TigrGffExon::initExonTypeNames()
{
  exonTypeNames["exon"]=ET_EXON;
  exonTypeNames["initial-exon"]=ET_INITIAL_EXON;
  exonTypeNames["final-exon"]=ET_FINAL_EXON;
  exonTypeNames["internal-exon"]=ET_INTERNAL_EXON;
  exonTypeNames["single-exon"]=ET_SINGLE_EXON;
}



void TigrGffExon::toGff(ostream &os)
{
  os << parent.getSubstrate() << "\t"
     << parent.getSource() << "\t"
     << exonType << "\t"
     << begin+1 << "\t"
     << end << "\t";
  if(hasScore) os << score; else os << ".";
  os << "\t" << parent.getStrand() << "\t";
  if(hasFrame) os << frame; else os << ".";
  os << "\ttransgrp=" << parent.getTranscriptId() << endl;
}



void TigrGffExon::changeExonType(ExonType e)
{
  exonType=e;
}


