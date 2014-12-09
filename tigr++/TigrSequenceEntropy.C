/****************************************************************
 TigrSequenceEntropy.C
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#include "TigrSequenceEntropy.H"
#include <iostream>
#include "TigrStringMap.H"
#include "TigrVector.H"
#include "DnaAlphabet.H"


double logOf2=log(2.0);
inline double lg(double x) {return log(x)/logOf2;}


inline int hashTableSize(int order)
{
  // these are all primes

  switch(order)
    {
    case 0: return 11; break;
    case 1: return 29; break;
    case 2: return 127; break;
    case 3: return 619; break;
    default:
    case 4: return 3121; break;
    }
}



double TigrSequenceEntropy::entropy(const Sequence &seq,int order,
				    double &maxEntropy)
{
  TigrString *str=seq.toString(DnaAlphabet::global);
  double H=entropy(*str,order,maxEntropy);
  delete str;
  return H;
}



double TigrSequenceEntropy::entropy(const TigrString &str,int order,
				    double &maxEntropy)
{
  int len=str.length();
  int gramSize=order+1;
  if(gramSize>=len) 
    throw TigrString("Order ")+order+" is too large for sequence of length "+len;
  int numWindows=len-gramSize+1;
  TigrStringMap<int> counts(hashTableSize(order));
  const char *p=str.c_str();
  int total=0;
  for(int i=0 ; i<numWindows ; ++i, ++p)
    {
      if(counts.isDefined(p,gramSize)) 
	++counts.lookup(p,gramSize);
      else 
	counts.lookup(p,gramSize)=1;
      ++total;
    }
  double dTotal=total, entropy=0;
  StringMapIterator<int> cur=counts.begin(), end=counts.end();
  for(; cur!=end ; ++cur)
    {
      int count=(*cur).second;
      double p=count/dTotal;
      //if(gramSize==1) cout<<TigrString((*cur).first,(*cur).len)<<"\t"<<count<<"/"<<total<<"="<<p<<endl;
      entropy-=p*lg(p);
    }
  maxEntropy=-lg(1.0/counts.size());
  if(entropy>maxEntropy) entropy=maxEntropy;
  return entropy;
}

