/****************************************************************
 TigrFastaWriter.H

 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_TigrFastaWriter_H
#define INCL_TigrFastaWriter_H
#include "TigrString.H"
#include <iostream>
#include "TigrRegex.H"
using namespace std;

class TigrFastaWriter
{
public:
  TigrFastaWriter(int width=60);
  void writeFasta(const TigrString &defline,const TigrString &sequence,
		  const TigrString &filename);
  void writeFastaFromCharPtr(const TigrString &defline,const char *sequence,
			 const TigrString &filename);
  void addToFasta(const TigrString &defline,const char *sequence,
		  ostream &);
  void addToFasta(const TigrString &defline,const TigrString &sequence,
		  ostream &);
private:
  int width;
  TigrRegex newline, greater;
};


#endif
