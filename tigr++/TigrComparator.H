/**************************************************************
TigrComparator.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrComparator_H
#define INCL_TigrComparator_H
using namespace std;

template <class T>
class TigrComparator
{
public:
  virtual bool equal(T &a,T &b)=0;
  virtual bool greater(T &a,T &b)=0;
  virtual bool less(T &a,T &b)=0;
};


template <class T>
class DirectComparator : public TigrComparator<T>
{
public:
  virtual bool equal(T &a,T &b) {return a==b;}
  virtual bool greater(T &a,T &b) {return a>b;}
  virtual bool less(T &a,T &b) {return a<b;}
};


#endif
