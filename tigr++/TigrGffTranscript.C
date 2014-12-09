#include <iostream>
#include "TigrVectorSorter.H"
#include "TigrGffTranscript.H"
using namespace std;



TigrGffTranscript::TigrGffTranscript(const TigrString &transcriptId,
			     const TigrString &substrate,
			     char strand,const TigrString &source)
  : transcriptId(transcriptId), strand(strand), begin(-1), end(-1),
    substrate(substrate), source(source), startCodon(NULL),
    stopCodon(NULL), hasScore(false)
{ 
  // ctor
}



TigrGffTranscript::~TigrGffTranscript()
{
  // dtor

  TigrVector<TigrGffExon*>::iterator cur=exons.begin(), end=exons.end();
  for(; cur!=end ; ++cur) delete *cur;
  delete startCodon;
  delete stopCodon;
}



ostream &operator<<(ostream &os,const TigrGffTranscript &transcript)
{
  transcript.printOn(os);
  return os;
}



TigrGffExon &TigrGffTranscript::getIthExon(int i)
{
  return *exons[i];
}



TigrGffFeature *TigrGffTranscript::getStartCodon()
{
  return startCodon;
}



TigrGffFeature *TigrGffTranscript::getStopCodon()
{
  return stopCodon;
}



char TigrGffTranscript::getStrand()
{
  return strand;
}



class ExonComparator : public TigrComparator<TigrGffExon*>
{
public:
  bool less(TigrGffExon *&a,TigrGffExon *&b) {return a->getBegin()<b->getBegin();}
  bool greater(TigrGffExon *&a,TigrGffExon *&b) {return a->getBegin()>b->getBegin();}
  bool equal(TigrGffExon *&a,TigrGffExon *&b) {return a->getBegin()==b->getBegin();}
};



const TigrString &TigrGffTranscript::getSource() const
{
  return source;
}



const TigrString &TigrGffTranscript::getSubstrate() const
{
  return substrate;
}



const TigrString &TigrGffTranscript::getTranscriptId() const
{
  return transcriptId;
}



double TigrGffTranscript::getScore() const
{
  return score;
}



int TigrGffTranscript::getBegin() const
{
  return begin;
}



int TigrGffTranscript::getEnd() const
{
  return end;
}



int TigrGffTranscript::getNumExons() const
{
  return exons.size();
}



void TigrGffTranscript:: printOn(ostream &os) const
{
  os << substrate << "\t"
     << source << "\t"
     << "transcript\t"
     << begin+1 << "\t"
     << end << "\t"
     << ".\t"
     << strand << "\t"
     << ".\ttransgrp=" << transcriptId
     << ";\tnumExons=" << getNumExons();
}



void TigrGffTranscript::addExon(TigrGffExon *exon)
{
  exons.push_back(exon);
  int exonBegin=exon->getBegin(), exonEnd=exon->getEnd();
  if(begin<0 || exonBegin<begin) begin=exonBegin;
  if(end<0 || exonEnd>end) end=exonEnd;
}



void TigrGffTranscript::setScore(double s)
{
  score=s;
  hasScore=true;
}



void TigrGffTranscript::setStartCodon(TigrGffFeature *s)
{
  startCodon=s;
}



void TigrGffTranscript::setStopCodon(TigrGffFeature *s)
{
  stopCodon=s;
}



void TigrGffTranscript::setStrand(char s)
{
  strand=s;
}



void TigrGffTranscript::sortExons()
{
  int numExons=exons.size();
  ExonComparator comp;
  TigrVectorSorter<TigrGffExon*> sorter(exons,comp);
  switch(strand)
    {
    case '+': 
      sorter.sortAscendInPlace(); 
      begin=exons[0]->getBegin();
      end=exons[numExons-1]->getEnd();
      break;
    case '-': 
      sorter.sortDescendInPlace(); 
      begin=exons[numExons-1]->getBegin();
      end=exons[0]->getEnd();
      break;
    default:  throw "strand is undefined in TigrGffTranscript::sortExons";
    }
}



void TigrGffTranscript::toGff(ostream &os)
{
  int n=exons.size();
  for(int i=0 ; i<n ; ++i)
    {
      TigrGffExon *exon=exons[i];
      exon->toGff(os);
    }

  if(startCodon) os << startCodon->toGff();
  if(stopCodon) os << stopCodon->toGff();

  os << substrate << "\t"
     << source << "\t"
     << "transcript\t"
     << begin+1 << "\t"
     << end << "\t";
  if(hasScore)
    os << score << "\t";
  else
    os << ".\t";
  os << strand << "\t"
     << ".\ttransgrp=" << transcriptId 
     << endl;
}



void TigrGffTranscript::setExonTypes()
{
  int numExons=exons.size();
  if(numExons==1)
    {
      exons[0]->changeExonType(ET_SINGLE_EXON);
      return;
    }

  exons[0]->changeExonType(ET_INITIAL_EXON);
  exons[numExons-1]->changeExonType(ET_FINAL_EXON);
  int nMinus1=numExons-1;
  for(int i=1 ; i<nMinus1 ; ++i)
    exons[i]->changeExonType(ET_INTERNAL_EXON);
}


