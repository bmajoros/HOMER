/**************************************************************
TigrVector.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrVector_H
#define INCL_TigrVector_H
using namespace std;
#include <vector>

template<class T>
class TigrVector : public vector<T>
{
public:
  TigrVector();
  TigrVector(int size);
  TigrVector<T> *getSubrange(int from,int to);
  void getSubrange(int from,int to,TigrVector<T> &);
  bool isMember(const T &);
  void append(const TigrVector<T> &);
  void push_front(const T &);
  void removeFirst(const T &);
  bool isEmpty() const {return empty();}
  T getFirst() {return *begin();}
  T getLast() {return (*this)[size()-1];}
  void resize(int n) {vector<T>::resize(n);}
  void cut(int index);
  void cut(int startIndex,int howMany);
  void insertByIndex(const T &,int beforeIndex);
  void reverse();
  inline void swap(int i,int j);
};

typedef TigrVector<int> IntVector;
typedef TigrVector<float> FloatVector;
typedef TigrVector<double> DoubleVector;
typedef TigrVector<char> CharVector;
typedef TigrVector<long> LongVector;


template<class T>
TigrVector<T>::TigrVector()
{
}



template<class T>
TigrVector<T>::TigrVector(int size) : vector<T>(size)
{
}



template<class T>
TigrVector<T> *TigrVector<T>::getSubrange(int from,int to)
{
  TigrVector<T> &newVector=*new TigrVector<T>;
  getSubrange(from,to,newVector);
  return &newVector;
}



template<class T>
bool TigrVector<T>::isMember(const T &x)
{
  typename TigrVector<T>::iterator cur=begin(), end=end();
  for(; cur!=end ; ++cur)
    if(*cur==x) 
      return true;
  return false;
}



template<class T>
void TigrVector<T>::append(const TigrVector<T> &other)
{
  typename TigrVector<T>::const_iterator cur=other.begin(), end=other.end();
  for( ; cur!=end ; ++cur)
    push_back(*cur);
}



template<class T>
void TigrVector<T>::getSubrange(int from,int to,TigrVector<T> &into)
{
  TigrVector<T> &self=*this;
  for(int i=from ; i<=to ; ++i)
    into.push_back(self[i]);
}



template<class T>
void TigrVector<T>::push_front(const T &x)
{
  insert(begin(),x);
}



template<class T>
void TigrVector<T>::removeFirst(const T &x)
{
  typename TigrVector<T>::iterator cur=begin(), end=this->end();
  for(; cur!=end ; ++cur)
    if(*cur==x) 
      {
	erase(cur);
	return;
      }
}



template<class T>
void TigrVector<T>::cut(int index)
{
  cut(index,1);
}



template<class T>
void TigrVector<T>::cut(int i,int howMany)
{
  TigrVector<T> &self=*this;
  int n=size();
  int src=i+howMany;
  for(; src<n ; ++i, ++src)
    self[i]=self[src];
  resize(n-howMany);
}



template<class T>
void TigrVector<T>::insertByIndex(const T &newT,int beforeIndex)
{
  int n=size();
  TigrVector<T> &self=*this;
  push_back(newT);
  for(int i=n ; i>beforeIndex ; --i) self[i]=self[i-1];
  self[beforeIndex]=newT;
}



template<class T>
void TigrVector<T>::swap(int i,int j)
{
  TigrVector<T> &self=*this;
  T t=self[i];
  self[i]=self[j];
  self[j]=t;
}



template<class T>
void TigrVector<T>::reverse()
{
  for(int i=0, j=size()-1 ; i<j ; ++i, --j)
    swap(i,j);
}


#endif
