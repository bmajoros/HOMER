/**************************************************************
TigrListQueue.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrListQueue_H
#define INCL_TigrListQueue_H
#include <list>
using namespace std;

template<class T>
class TigrListQueue
{
public:
  TigrListQueue();
  T &back();
  T &front();
  typename list<T>::iterator begin();
  typename list<T>::iterator end();
  typename list<T>::reverse_iterator rbegin(); 
  typename list<T>::reverse_iterator rend(); 
  void enqueue(const T &);
  T dequeue() {T &t=back();pop_back();return t;}
  void pop_back();
  void pop_front();
  int size();
  void clear();
  bool isEmpty() {return theList.empty();}
private:
  list<T> theList;
};



template<class T>
int TigrListQueue<T>::size()
{
  return theList.size();
}



template<class T>
void TigrListQueue<T>::clear()
{
  theList.clear();
}



template<class T>
T &TigrListQueue<T>::back()
{
  return theList.back();
}



template<class T>
T &TigrListQueue<T>::front()
{
  return theList.front();
}



template<class T>
TigrListQueue<T>::TigrListQueue() 
{
   
}



template<class T>
typename list<T>::iterator TigrListQueue<T>::begin()
{
  return theList.begin();
}



template<class T>
typename list<T>::iterator TigrListQueue<T>::end()
{
  return theList.end();
}



template<class T>
typename list<T>::reverse_iterator TigrListQueue<T>::rbegin()
{
  return theList.rbegin();
}



template<class T>
typename list<T>::reverse_iterator TigrListQueue<T>::rend()
{
  return theList.rend();
}



template<class T>
void TigrListQueue<T>::enqueue(const T &t)
{
  theList.push_back(t);
}



template<class T>
void TigrListQueue<T>::pop_back()
{
  theList.pop_back();
}



template<class T>
void TigrListQueue<T>::pop_front()
{
  theList.pop_front();
}

#endif
