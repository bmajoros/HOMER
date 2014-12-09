/**************************************************************
TigrVectorPrinter.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrVectorPrinter_H
#define INCL_TigrVectorPrinter_H
#include <iostream>
#include "TigrVector.H"
using namespace std;

template<class T>
ostream &operator<<(ostream &os,TigrVector<T> &v)
{
  os << '[';
  int i, n=v.size();
  for(i=0 ; i<n ; ++i)
    os << v[i] << ',';
  os << ']';
}

template<class T>
void printVectorUsingIndex(ostream &os,TigrVector<T> &v,TigrVector<int> &I)
{
  os << '[';
  int i, n=v.size();
  for(i=0 ; i<n ; ++i)
    os << v[I[i]] << ',';
  os << ']';
 
}

#endif
