#include <iostream>
#include "TigrVectorSorter.H"
#include "TigrGffReader.H"
using namespace std;


bool TranscriptComparator::less(TigrGffTranscript *&a,TigrGffTranscript *&b) 
{
  return a->getBegin()<b->getBegin();
}



bool TranscriptComparator::greater(TigrGffTranscript *&a,
				   TigrGffTranscript *&b) 
{
  return a->getBegin()>b->getBegin();
}



bool TranscriptComparator::equal(TigrGffTranscript *&a,TigrGffTranscript *&b) 
{
  return a->getBegin()==b->getBegin();
}



TigrGffReader::TigrGffReader(const TigrString &filename)
  : commentPattern("\\s*#.*"), transgrpRegex("transgrp=([^;]+)")
{
  if(!file.open(filename.c_str(),"rb"))
    throw TigrString("Can't open file ")+filename;

  exonTypes.insert("exon");
  exonTypes.insert("initial-exon");
  exonTypes.insert("internal-exon");
  exonTypes.insert("final-exon");
  exonTypes.insert("single-exon");
}



TigrGffFeature *TigrGffReader::nextFeature()
{
  while(true)
    {
      TigrString line=file.readLine();
      if(file.eof()) break;
      if(line.length()==0) continue;
      if(commentPattern.match(line)) continue;
      return new TigrGffFeature(line);
    }
  return NULL;
}



TigrVector<TigrGffTranscript*> *TigrGffReader::loadTranscripts()
{
  // Read the features fromthe GFF file
  TigrMap<TigrString,TigrGffTranscript*> transHash;
  while(TigrGffFeature *f=nextFeature())
    {
      // If this feature is an exon...
      if(f->hasExtraFields() && exonTypes.isMember(f->getFeatureType()))
	{
	  // Parse out transcript ID
	  TigrString transcriptId=f->getExtraFields()[0];
	  if(transgrpRegex.search(transcriptId)) 
	    transcriptId=transgrpRegex[1];
	  if(!transHash.isDefined(transcriptId))
	    transHash[transcriptId]=
	      new TigrGffTranscript(transcriptId,f->getSubstrate(),
				f->getStrand(),f->getSource());

	  // Add this exon to the appropriate transcript
	  TigrGffTranscript *transcript=transHash[transcriptId];
	  TigrGffExon *exon=new TigrGffExon(*f,*transcript);
	  transcript->addExon(exon);
	}
      delete f;
    }

  // Sort exons and resolve exon types; also convert the hash table
  // to a vector of transcripts
  TigrVector<TigrGffTranscript*> *transcriptList=
    new TigrVector<TigrGffTranscript*>;
  TigrMap<TigrString,TigrGffTranscript*>::iterator cur=transHash.begin(),
    end=transHash.end();
  for(; cur!=end ; ++cur)
    {
      TigrGffTranscript *transcript=(*cur).second;
      transcript->sortExons();
      transcript->setExonTypes();
      transcriptList->push_back(transcript);
    }

  // Sort the transcripts by position
  TranscriptComparator comp;
  TigrVectorSorter<TigrGffTranscript*> sorter(*transcriptList,comp);
  sorter.sortAscendInPlace();

  return transcriptList;
}



TigrMap<TigrString,TigrGffReader::TranscriptList*> *
  TigrGffReader::loadByContig()
{
  TigrMap<TigrString,TranscriptList*> &M=*
    new TigrMap<TigrString,TranscriptList*>;
  TranscriptList *allTranscripts=loadTranscripts();
  TigrGffReader::TranscriptList::iterator cur=allTranscripts->begin(),
    end=allTranscripts->end();
  for(; cur!=end ; ++cur)
    {
      TigrGffTranscript *transcript=*cur;
      if(!transcript) break;
      const TigrString &substrate=transcript->getSubstrate();
      if(!M.isDefined(substrate)) M[substrate]=new TranscriptList;
      TranscriptList *transcripts=M[substrate];
      transcripts->push_back(transcript);
    }
  delete allTranscripts;
  return &M;
}



