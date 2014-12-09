/**************************************************************
TigrSet.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrSet_H
#define INCL_TigrSet_H
#include <set>
#include <algorithm>
using namespace std;

template< class T, class C = less<T> >
class TigrSet : public set<T,C>
{
public:
  TigrSet() {}
  bool operator+=(const TigrSet<T,C> &);
  bool operator+=(const T &);
  bool isMember(const T &);
  void remove(const T &t);
  TigrSet<T,C> *intersectWith(TigrSet<T,C> &);
  void intersectWith(TigrSet<T,C> &,TigrSet<T,C> &destination);
  inline bool isEmpty() const {return empty();}
};



template<class T, class C>
void TigrSet<T,C>::intersectWith(TigrSet<T,C> &other,
				 TigrSet<T,C> &destination)
{
  set_intersection(begin(),end(),other.begin(),other.end(),
		   inserter(destination,destination.begin()));
}



template<class T, class C>
TigrSet<T,C> *TigrSet<T,C>::intersectWith(TigrSet<T,C> &other)
{
  TigrSet<T,C> *newSet=new TigrSet<T,C>;
  intersectWith(other,*newSet);
  return newSet;
}



template<class T, class C>
bool TigrSet<T,C>::isMember(const T &t)
{
  return find(t)!=end();
}



template<class T, class C>
bool TigrSet<T,C>::operator+=(const T &other)
{
  if(!isMember(other))
    {
      insert(other);
      return true;
    }
  return false;
}



template<class T, class C>
bool TigrSet<T,C>::operator+=(const TigrSet<T,C> &other)
{
  bool changes=false;
  typename TigrSet<T,C>::iterator cur=other.begin(), end=other.end();
  for(; cur!=end ; ++cur)
    if(!isMember(*cur))
      {
	insert(*cur);
	changes=true;
      }
  return false;
}



template<class T, class C>
void TigrSet<T,C>::remove(const T &t) 
{
  erase(t);
}

#endif
