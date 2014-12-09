/****************************************************************
 TigrStack.H

 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_TigrStack_H
#define INCL_TigrStack_H
#include "TigrVector.H"
using namespace std;

template<class T>
class TigrStack
{
public:
  TigrStack();
  TigrStack(TigrVector<T> &);
  virtual ~TigrStack();
  void push(const T &);
  T pop(); // removes & returns the top element!
  void discardTop(); // removes top element
  const T &top();
  int size() const;
  bool isEmpty() const;
  typedef typename TigrVector<T>::iterator iterator;
  typedef typename TigrVector<T>::reverse_iterator reverse_iterator;
  iterator begin();
  iterator end();
  reverse_iterator rbegin();
  reverse_iterator rend();
  TigrVector<T> &asVector() {return theStack;}
private:
  TigrVector<T> &theStack;
  bool shouldDelete;
};


template<class T>
TigrStack<T>::TigrStack()
  : theStack(*new TigrVector<T>),
    shouldDelete(true)
{
}



template<class T>
TigrStack<T>::TigrStack(TigrVector<T> &v)
  : theStack(v),
    shouldDelete(false)
{
}



template<class T>
TigrStack<T>::~TigrStack()
{
  if(shouldDelete) delete &theStack;
}



template<class T>
void TigrStack<T>::push(const T &t)
{
  theStack.push_back(t);
}



template<class T>
T TigrStack<T>::pop()
{
  T t=theStack.back();
  theStack.resize(theStack.size()-1);
  return t;
}



template<class T>
void TigrStack<T>::discardTop()
{
  theStack.resize(theStack.size()-1);
}



template<class T>
const T &TigrStack<T>::top()
{
  return theStack.back();
}



template<class T>
int TigrStack<T>::size() const
{
  return theStack.size();
}



template<class T>
bool TigrStack<T>::isEmpty() const
{
  return theStack.size()==0;
}



template<class T>
typename TigrStack<T>::iterator TigrStack<T>::begin()
{
  return theStack.begin();
}



template<class T>
typename TigrStack<T>::iterator TigrStack<T>::end()
{
  return theStack.end();
}



template<class T>
typename TigrStack<T>::reverse_iterator TigrStack<T>::rbegin()
{
  return theStack.rbegin();
}



template<class T>
typename TigrStack<T>::reverse_iterator TigrStack<T>::rend()
{
  return theStack.rend();
}


#endif
