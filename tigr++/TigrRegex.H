/**************************************************************
TigrRegex.H
bmajoros@tigr.org 1/1/2003

This code utilizes the GNU regular expression library, which is
copyrighted property of the Free Software Foundation, Inc., and
is distributed under the GNU General Public License.

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrRegex_H
#define INCL_TigrRegex_H
#include "TigrString.H"
using namespace std;

extern "C" {
#include "regex.h"
}

class TigrRegex
{
public:
  TigrRegex(const TigrString &regex);
  ~TigrRegex();
  TigrString operator[](int whichSubexpression);  
  TigrString getEntireMatch();
  TigrString substitute(const TigrString &substrate,
			const TigrString &replacement);
  TigrVector<TigrString*> *split(const TigrString &substrate);
  bool match(const TigrString &);
  bool search(const TigrString &);
  int getNumSubexpressions() const;
  static TigrString substitute(const TigrString &regex,
			       const TigrString &replacement,
			       const TigrString &substrate);
  static TigrVector<TigrString*> *split(const TigrString &regex,
				    const TigrString &substrate);
  static bool match(const TigrString &regex,const TigrString &substrate);
  static bool search(const TigrString &regex,const TigrString &substrate);
  void getMatchIndices(int regNum,int &start,int &end);
private:
  re_registers reg;
  re_pattern_buffer patbuf;
  int numOfSubexpress;
  TigrString regex;
  const TigrString *substrate;
  TigrString preprocess(const TigrString &);
  void compile();
};

#endif
