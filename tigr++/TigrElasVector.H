/**************************************************************
TigrElasVector.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrElasVector_H
#define INCL_TigrElasVector_H
#include "TigrVector.H"
using namespace std;

template<class T>
class TigrElasVector
{
public:
  TigrElasVector(const T &nullInitializer);
  T &operator[](int);
  int size();
private:
  const T &nullInitializer;
  TigrVector<T> theVector;
};

class TigrElasDblVect : public TigrElasVector<double>
{
public:
  TigrElasDblVect();
};

class TigrElasFloatVect : public TigrElasVector<float>
{
public:
  TigrElasFloatVect();
};

class TigrElasIntVect : public TigrElasVector<int>
{
public:
  TigrElasIntVect();
};

class TigrElasCharVect : public TigrElasVector<char>
{
public:
  TigrElasCharVect();
};

class TigrElasLongVect : public TigrElasVector<long>
{
public:
  TigrElasLongVect();
};



TigrElasDblVect() : TigrElasVector<double>(0) 
{
}



TigrElasFloatVect() : TigrElasVector<float>(0) 
{
}



TigrElasIntVect() : TigrElasVector<int>(0) 
{
}



TigrElasCharVect() : TigrElasVector<char>('\0') 
{
}



TigrElasLongVect() : TigrElasVector<long>(0) 
{
}



template<class T>
T &TigrElasVector<T>::operator[](int index)
{
  int sz=theVector.size();
  if(index>=sz)
    for(int i=sz ; i<=index ; ++i)
      theVector.push_back(nullInitializer);
  return theVector[index];
}



template<class T>
TigrElasVector<T>::TigrElasVector(const T &nullInitializer)
  : nullInitializer(nullInitializer)
{
}



template<class T>
int TigrElasVector<T>::size()
{
  return theVector.size();
}

#endif
