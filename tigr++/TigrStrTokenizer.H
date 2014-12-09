/**************************************************************
TigrStrTokenizer.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrStrTokenizer_H
#define INCL_TigrStrTokenizer_H
#include "TigrVector.H"
#include "TigrString.H"
using namespace std;

#define DEFAULT_WHITESPACE " \t\n\r"

typedef pair< TigrVector<TigrString>, TigrVector<TigrString> > StrVectPair;

class TigrStrTokenizer 
{
public:
  TigrStrTokenizer(const char *wspace=DEFAULT_WHITESPACE);
  TigrStrTokenizer(const char *source,
		  const char *wspace=DEFAULT_WHITESPACE);
  TigrStrTokenizer(const TigrString &source,
		  const char *wspace=DEFAULT_WHITESPACE);
  TigrVector<TigrString*> *getAllTokens();
  TigrVector<TigrString> *getTokenStrings();
  bool hasMoreTokens();
  const char *nextToken();
  StrVectPair *getTokensAndSeparators(const TigrString &);
  static TigrVector<TigrString*> *tokenize(const char *source,
		  const char *wspace=DEFAULT_WHITESPACE);
  void tokenizeThis(const TigrString &s);
  void tokenizeThis(const char *p);
private:
  enum WhereType 
    { 
      IN_TOKEN, 
      IN_SPACE, 
      AT_START 
    };
  const char *p;
  const char *source;
  const char *wspace;
  TigrString buffer;
  bool isWhiteSpace(char c);
  void skipWhiteSpace();
};

#endif
