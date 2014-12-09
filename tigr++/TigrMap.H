/**************************************************************
TigrMap.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrMap_H
#define INCL_TigrMap_H
#include <map>
using namespace std;

template<class T>
struct DefaultMapComparator
{
  bool operator()(const T &a,const T &b){return a<b;}
};

template<class KeyType,class ValueType,
	 class Comp=DefaultMapComparator<KeyType> >
class TigrMap : public map<KeyType,ValueType,Comp>
{
public:
  bool isDefined(const KeyType &key);
  void undefine(const KeyType &key) {erase(key);}
  inline int numElements();
  inline bool isEmpty() {return numElements()==0;}
};



template<class KeyType,class ValueType,class Comp>
bool TigrMap<KeyType,ValueType,Comp>::isDefined(const KeyType &key)
{
  return find(key)!=end();
}



template<class KeyType,class ValueType,class Comp>
inline int TigrMap<KeyType,ValueType,Comp>::numElements() 
{
  return size();
}


#endif
