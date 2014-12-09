/****************************************************************
 TigrHugeFastaReader.C

 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/



#include "TigrHugeFastaReader.H"
#include <iostream>
using namespace std;


TigrHugeFastaReader::TigrHugeFastaReader(const TigrString &filename,
					 int chunkSize,
					 int overlapSize)
  : file(filename), chunkSize(chunkSize),
    overlapSize(overlapSize), buffer(new char[chunkSize+1])
{
  if(overlapSize<0 || overlapSize>=chunkSize)
    throw "Bad overlap size in TigrHugeFastaReader";

  remainderSize=chunkSize-overlapSize;
  overlap=buffer+remainderSize;
  remainder=buffer+overlapSize;
  end=buffer+chunkSize;
  defline=file.getline();
  loadFirstChunk();
}



TigrHugeFastaReader::~TigrHugeFastaReader()
{
  delete [] buffer;
}



const char *TigrHugeFastaReader::getBuffer()
{
  return buffer;
}



void TigrHugeFastaReader::loadFirstChunk()
{
  int bytesRead=file.read(chunkSize,buffer);
  buffer[bytesRead]='\0';
  filterSpaces(buffer);
}



bool TigrHugeFastaReader::loadNextChunk()
{
  // First, move overlap region to beginning of buffer
  memmove(buffer,overlap,overlapSize);
  
  // Now load remainder of chunk into the rest of the buffer
  int bytesRead=file.read(remainderSize,remainder);
  remainder[bytesRead]='\0';
  filterSpaces(remainder);
  return !file.eof();
}



bool TigrHugeFastaReader::isMore()
{
  return !file.eof();
}



void TigrHugeFastaReader::filterSpaces(char *buf)
{
  // Compact the sequence to remove any spaces
  char *src=buf, *dest=buf;
  for( ; *src ; ++src)
    {
      char c=*src;
      if(isspace(c)) continue;
      if(src>dest) *dest=c;
      ++dest;
    }

  // Now read in any additional bases to fill up the rest of the buffer
  while(dest!=end && !file.eof())
    {
      if(!file.read(1,dest)) break;
      if(!isspace(*dest)) ++dest;
    }
  *dest='\0';
}



const TigrString &TigrHugeFastaReader::getDefline()
{
  return defline;
}


