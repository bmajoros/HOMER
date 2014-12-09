/****************************************************************
 Histogram.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_Histogram_H
#define INCL_Histogram_H
#include "TigrArray1D.H"
#include <iostream>
#include <fstream>

template<class T>
class Histogram
{
  TigrArray1D<T> bins;
  T minValue, maxValue, binSize;

  void load(const TigrString &filename);
public:
  typedef double PROBABILITY;
  Histogram(T minValue,T maxValue,int numBins,float pseudocount);
  Histogram(const TigrString &filename);
  void addCounts(TigrVector<T> &);
  void normalize();
  PROBABILITY getP(T) const;
  PROBABILITY getBin(int);
  int getNumBins() const;
  T getMin() const;
  T getMax() const;
  void getBinBoundaries(int binNum,T &binMin,T &binMax) const;
  void divideBy(const Histogram<T> &other);
  void save(const TigrString &filename);
  void useLogs();
  T sum() const;
};



template<class T>
Histogram<T>::Histogram(T minValue,T maxValue,int numBins,float pseudocount)
  : bins(numBins), 
    minValue(minValue),
    maxValue(maxValue),
    binSize((maxValue-minValue)/numBins)
{
  // ctor

  bins.setAllTo(pseudocount);
}



template<class T>
Histogram<T>::Histogram(const TigrString &filename)
  : bins(0)
{
  // ctor

  load(filename);
}



template<class T>
void Histogram<T>::addCounts(TigrVector<T> &v)
{
  typename TigrVector<T>::iterator cur=v.begin(), end=v.end();
  for(; cur!=end ; ++cur)
    {
      const T &t=*cur;
      int binNum=int((t-minValue)/binSize);
      if(binNum<0) binNum=0;
      else if(binNum>=bins.size()) binNum=bins.size()-1;
      bins[binNum]+=(T)1;
    }
}



template<class T>
void Histogram<T>::normalize()
{
  T total=0;
  int numBins=bins.size();
  for(int i=0 ; i<numBins ; ++i)
    total+=bins[i];
  for(int i=0 ; i<numBins ; ++i)
    bins[i]/=total;
}



template<class T>
double Histogram<T>::getP(T t) const
{
  int binNum=int((t-minValue)/binSize);
  if(binNum<0) binNum=0;
  else if(binNum>=bins.size()) binNum=bins.size()-1;
  return bins[binNum];
}



template<class T>
int Histogram<T>::getNumBins() const
{
  return numBins;
}



template<class T>
T Histogram<T>::getMin() const
{
  return minValue;
}



template<class T>
T Histogram<T>::getMax() const
{
  return maxValue;
}



template<class T>
void Histogram<T>::getBinBoundaries(int binNum,T &binMin,T &binMax) const
{
  binMin=binNum*binSize;
  binMax=binMin+binSize;
}



template<class T>
double Histogram<T>::getBin(int i)
{
  return bins[i];
}



template<class T>
void Histogram<T>::divideBy(const Histogram<T> &other)
{
  int n=bins.size();
  for(int i=0 ; i<n ; ++i)
    bins[i]/=other.bins[i];
}



template<class T>
void Histogram<T>::save(const TigrString &filename)
{
  ofstream os(filename.c_str());
  if(!os.good()) throw TigrString("Can't write into file: ")+filename;
  int n=bins.size();
  os<<minValue<<"\t"<<maxValue<<"\t"<<n<<"\t"<<binSize<<endl;
  T binBegin=minValue;
  for(int i=0 ; i<n ; ++i)
    {
      T binEnd=binBegin+binSize;
      os<<binBegin<<"\t"<<binEnd<<"\t"<<bins[i]<<endl;
      binBegin=binEnd;
    }
}



template<class T>
void Histogram<T>::load(const TigrString &filename)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw TigrString("Can't open file: ")+filename;
  int numBins;
  is>>minValue>>maxValue>>numBins>>binSize;
  bins.resize(numBins);
  T begin, end;
  for(int i=0 ; i<numBins ; ++i)
    is>>begin>>end>>bins[i];
}



template<class T>
void Histogram<T>::useLogs()
{
  int n=bins.size();
  for(int i=0 ; i<n ; ++i)
    bins[i]=log(bins[i]);
}



template<class T>
T Histogram<T>::sum() const
{
  T sum=0;
  int n=bins.size();
  for(int i=0 ; i<n ; ++i)
    sum+=bins[i];
  return sum;
}

#endif
