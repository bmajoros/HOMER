/****************************************************************
 TigrFastaWriter.C

 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#include "TigrFastaWriter.H"
#include <fstream>
using namespace std;


TigrFastaWriter::TigrFastaWriter(int width)
  : width(width), newline("\n$"), greater("^\\s*>")
{
}



void TigrFastaWriter::writeFasta(const TigrString &defline,
				 const TigrString &sequence,
				 const TigrString &filename)
{
  ofstream os(filename.c_str());
  addToFasta(defline,sequence.c_str(),os);
}



void TigrFastaWriter::writeFastaFromCharPtr(const TigrString &defline,
					    const char *sequence,
					    const TigrString &filename)
{
  ofstream os(filename.c_str());
  addToFasta(defline,sequence,os);
}



void TigrFastaWriter::addToFasta(const TigrString &defline,
				 const TigrString &sequence,
				 ostream &os)
{
  addToFasta(defline,sequence.c_str(),os);
}



void TigrFastaWriter::addToFasta(const TigrString &def,
				 const char *sequence,
				 ostream &os)
{
  TigrString defline=def;
  if(newline.search(defline)) defline.chop();
  if(!greater.search(defline)) defline=TigrString(">")+defline;

  os << defline << endl;
  int length=strlen(sequence);
  int numLines=int(length/width);
  if(length%width) ++numLines;
  int start=0;

  char *line=new char[width+1];
  line[width]='\0';
  for(int i=0 ; i<numLines ; ++i)
    {
      strncpy(line,sequence+start,width);
      os << line << endl;
      start+=width;
    }
  delete [] line;
  if(length==0) os << endl;
}


