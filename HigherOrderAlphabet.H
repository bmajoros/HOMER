/****************************************************************
 HigherOrderAlphabet.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_HigherOrderAlphabet_H
#define INCL_HigherOrderAlphabet_H

#include "tigr++/Alphabet.H"
#include "tigr++/TigrMap.H"
#include "tigr++/TigrString.H"
#include "tigr++/TigrArray2D.H"

/******************************************************************
  HigherOrderAlphabet provides a mapping between ngrams and unique
  integer identifiers, so that the "symbols" in the higher-order
  alphabet are actually ngrams.  This class can also be used to
  initialize an "ngram array", which is a 2-dimensional array in
  which the first dimension is as long as the sequence and the
  second dimension is interpreted as the order.  Thus, given a
  sequence S, once we have initialized the ngram-array A for an
  Nth order HMM, the higher-order symbol at location A[i][j] is
  the integer identifier of the ngram ending at position i in
  the sequence and having order j (or length j+1).
******************************************************************/
class HigherOrderAlphabet
{
  Alphabet &alphabet;
  int order, numLabels;
  TigrMap<TigrString,int> symbolMap;
  TigrVector<int> firstOfOrder;

  void initSymbolMap();
public:
  HigherOrderAlphabet(Alphabet &,int order);
  void initNgramArray(TigrArray2D<int> &ngramArray, // index==order
		      const TigrString &seqStr); 
  int lookup(const TigrString &);
  int getOrder() const {return order;}
  int size();
  int firstSymbolOfOrder(int order);
  int lastSymbolOfOrder(int order);
  int getOrder(int ngram);
};


#endif
