/****************************************************************
 FixedSizePriorityList.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_FixedSizePriorityList_H
#define INCL_FixedSizePriorityList_H
using namespace std;
#include "TigrList.H"
#include "TigrComparator.H"
#include "TigrIterator.H"


template <class T>
class FSPLIterator : public TigrIterator<T>
{
  typename TigrList<T>::iterator iter;
public:
  FSPLIterator(typename TigrList<T>::iterator);
  FSPLIterator() {}
  virtual void operator++();
  virtual void operator++(int);
  virtual T &operator*();
  virtual bool operator==(const TigrIterator<T> &) const;
  virtual bool operator!=(const TigrIterator<T> &i) const
    {return !(i==*this);}
  virtual TigrIterator<T> &clone();
  operator typename TigrList<T>::iterator() {return iter;}
};



template <class T>
class FixedSizePriorityList
{
  TigrList<T> Q;
  TigrComparator<T> &comp;
  int capacity;
  void force_insert(T);
public:
  FixedSizePriorityList(int capacity,
			TigrComparator<T>& = *new DirectComparator<T>());
  int getCapacity() const;
  int size() const;
  bool empty() const;
  bool insert(T); // <---can cause memory leak if you're not careful! ###
  bool push_back(T t) {return insert(t);}
  T &getMin();
  T &getMax();
  typedef FSPLIterator<T> iterator;
  iterator begin();
  iterator end();
  void erase(iterator);
  void clear() {Q.clear();}
  TigrList<T> &getList() {return Q;}
};



template <class T>
FixedSizePriorityList<T>::FixedSizePriorityList(int capacity,
						TigrComparator<T> &comp)
  : capacity(capacity),
    comp(comp)
{
}



template <class T>
int FixedSizePriorityList<T>::getCapacity() const
{
  return capacity;
}



template <class T>
int FixedSizePriorityList<T>::size() const
{
  return Q.size();
}



template <class T>
bool FixedSizePriorityList<T>::empty() const
{
  return Q.size()==0;
}



template <class T>
bool FixedSizePriorityList<T>::insert(T t)
{
  // If the list is not yet full, just insert the new element:
  if(Q.size()<capacity)
    {
      force_insert(t);
      return true;
    }

  // The list is at capacity.  If the new element is smaller than
  // even the minimum element, don't insert it:
  if(comp.less(t,getMin())) return false;

  // The new element deserves insertion.  Insert it and then delete
  // the smallest element from the list.
  force_insert(t);
  Q.erase(Q.begin());
  return true;
}



template <class T>
void FixedSizePriorityList<T>::force_insert(T t)
{
  bool inserted=false;
  typename TigrList<T>::iterator cur=Q.begin(), end=Q.end();
  for(; cur!=end ; ++cur)
    if(comp.less(t,*cur))
      {
	Q.insert(cur,t);
	inserted=true;
	break;
      }
  if(!inserted) Q.push_back(t);
}



template <class T>
T &FixedSizePriorityList<T>::getMin()
{
  return *Q.begin();
}



template <class T>
T &FixedSizePriorityList<T>::getMax()
{
  return *Q.rbegin();
}



template <class T>
typename FixedSizePriorityList<T>::iterator FixedSizePriorityList<T>::begin()
{
  return FSPLIterator<T>(Q.begin());
}



template <class T>
typename FixedSizePriorityList<T>::iterator FixedSizePriorityList<T>::end()
{
  return FSPLIterator<T>(Q.end());
}



template <class T>
void FixedSizePriorityList<T>::erase(iterator i)
{
  Q.erase(i);
}



//================================================================
//                   FSPLIterator methods
//================================================================

template <class T>
FSPLIterator<T>::FSPLIterator(typename TigrList<T>::iterator iter)
  : iter(iter)
{
}



template <class T>
void FSPLIterator<T>::operator++()
{
  ++iter;
}



template <class T>
void FSPLIterator<T>::operator++(int)
{
  iter++;
}



template <class T>
T &FSPLIterator<T>::operator*()
{
  return *iter;
}



template <class T>
bool FSPLIterator<T>::operator==(const TigrIterator<T> &other) const
{
  return iter==dynamic_cast<const FSPLIterator<T>&>(other).iter;
}



template <class T>
TigrIterator<T> &FSPLIterator<T>::clone()
{
  return *new FSPLIterator(iter);
}


#endif
