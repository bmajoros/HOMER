/****************************************************************
 HigherOrderAlphabet.C
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#include "HigherOrderAlphabet.H"
#include <iostream>
#include "tigr++/NthOrderStringIterator.H"
using namespace std;


HigherOrderAlphabet::HigherOrderAlphabet(Alphabet &alphabet,int order)
  : alphabet(alphabet),
    order(order)
{
  initSymbolMap();
}



void HigherOrderAlphabet::initSymbolMap()
{
  int symbolId=0;
  for(int i=1 ; i<=order+1 ; ++i)
    {
      firstOfOrder.push_back(symbolId);
      NthOrderStringIterator iter(i,alphabet);
      while(!iter.done())
	{
	  TigrString label=iter.getNextString();
	  symbolMap[label]=symbolId;
	  ++symbolId;
	}
    }
  numLabels=symbolId;
  firstOfOrder.push_back(numLabels);
}



void HigherOrderAlphabet::initNgramArray(TigrArray2D<int> &ngramArray,
					 const TigrString &seqStr)
{
  int seqLen=seqStr.length();
  for(int pos=0 ; pos<seqLen ; ++pos)
    {
      int begin=pos-order;
      if(begin<0) begin=0;
      int ngramLen=pos-begin+1;
      TigrString fullNgram=seqStr.substr(begin,ngramLen);
      int symbolId=symbolMap[fullNgram];
      for(int i=ngramLen-1 ; i<=order ; ++i)
	ngramArray[pos][i]=symbolId;
      for(int i=1 ; i<ngramLen ; ++i)
	ngramArray[pos][ngramLen-i-1]=
	  symbolMap[fullNgram.substr(i,ngramLen-i)];
    }
}



int HigherOrderAlphabet::size()
{
  return numLabels;
}



int HigherOrderAlphabet::lookup(const TigrString &s)
{
  if(symbolMap.isDefined(s)) return symbolMap[s];
  throw s+" is not defined in HigherOrderAlphabet";
}




int HigherOrderAlphabet::firstSymbolOfOrder(int order)
{
  return firstOfOrder[order];
}



int HigherOrderAlphabet::lastSymbolOfOrder(int order)
{
  return firstOfOrder[order+1]-1;
}



int HigherOrderAlphabet::getOrder(int ngram)
{
  int n=firstOfOrder.size();
  for(int i=0 ; i<n ; ++i)
    if(ngram<firstOfOrder[i])
      return i-1;
  return n-1;
}


