/****************************************************************
 TigrHugeFastaReader.H

 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/

#ifndef INCL_TigrHugeFastaReader_H
#define INCL_TigrHugeFastaReader_H
using namespace std;

#include "TigrString.H"
#include "TigrFile.H"

// THIS CLASS HAS BEEN FULLY VERIFIED -- IT IS ALMOST *CERTAINLY*
// BUG FREE!

/****************************************************************
  This class can be used to read arbitrarily large FASTA files.
  The file is read in chunks.  You specify a chunk size (C) and an
  overlap size (O), and this class gives you successive chunks in
  the file.  The first chunk will be of the full chunk size C and
  is accessible via getBuffer(), but all subsequent chunks will be
  of size C-O and are accessed through getRemainder().  However,
  for the default overlap size of zero, you can just use getBuffer()
  for all chunks.  Incidentally, the buffer pointer that is returned
  does not change after calling loadNextChunk(), so you can just
  call it once and store it locally; loadNextChunk() will change
  the contents of that buffer for you each time.  The same is true
  of getRemainder(); it is just an offset into the buffer.

  The overlap feature is useful if you want to slide a window over
  the sequence without treating the ends of the chunks differently.
 */
class TigrHugeFastaReader
{
public:
  TigrHugeFastaReader(const TigrString &filename,int chunkSize,
		      int overlapSize=0);
  virtual ~TigrHugeFastaReader();
  const char *getBuffer();
  const char *getRemainder() {return remainder;}
  bool loadNextChunk();
  bool isMore();
  const TigrString &getDefline();
private:
  TigrString defline;
  TigrFile file;
  char *buffer, *overlap, *remainder, *end;
  int chunkSize, overlapSize, remainderSize;
  void loadFirstChunk();
  void filterSpaces(char *);
};


#endif
