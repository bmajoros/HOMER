/**************************************************************
TigrWindowedQueue.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrWindowedQueue_H
#define INCL_TigrWindowedQueue_H
#include <list>
using namespace std;

template<class T>
class TigrWindowedQueue
{
public:
  TigrWindowedQueue();
  bool isEmpty();
  list<T> &getReservoir();
  list<T> &getWindow();
  void enqueue(T);
  void shift();
private:
  list<T> reservoir, window;
};



template<class T>
TigrWindowedQueue<T>::TigrWindowedQueue()
{
   
}



template<class T>
bool TigrWindowedQueue<T>::isEmpty()
{
  return window.empty() && reservoir.empty();
}



template<class T>
list<T> &TigrWindowedQueue<T>::getReservoir()
{
  return reservoir;
}



template<class T>
list<T> &TigrWindowedQueue<T>::getWindow()
{
  return window;
}



template<class T>
void TigrWindowedQueue<T>::enqueue(T t)
{
  window.push_back(t);
}



template<class T>
void TigrWindowedQueue<T>::shift()
{
  reservoir.push_back(window.front());
  window.pop_front();
}

#endif
