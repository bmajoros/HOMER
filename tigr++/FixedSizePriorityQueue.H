/****************************************************************
 FixedSizePriorityQueue.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/

#ifndef INCL_FixedSizePriorityQueue_H
#define INCL_FixedSizePriorityQueue_H
using namespace std;
#include "TigrFastBinTree.H"

/*============================================================
  A priority queue which can hold a maximum of N elements.  If
  an element is inserted into a full queue, either the minimum
  element in the queue will be deleted (if the new element is
  larger), or the new element will be rejected.  Based on a
  TigrFastBinTree (fast binary tree in which each node has a
  pointer to its parent for fast deletion).

  Class T must define an operator<() and an operator==().
  ============================================================*/

template <class T>
class FixedSizePriorityQueue
{
public:
  FixedSizePriorityQueue(int capacity,
			 TigrComparator<T>& = *new DirectComparator<T>());
  int getCapacity() const;
  int getNumElements() const;
  bool isEmpty() const;
  bool insert(T);
  T &peekMin() const;
  T &peekMax() const;
  typedef TigrFastBinTree<T>::iterator iterator;
  iterator begin();
  iterator end();
private:
  int capacity;
  TigrFastBinTree<T> tree;
};



template <class T>
FixedSizePriorityQueue<T>::FixedSizePriorityQueue(int capacity,
						  TigrComparator<T> &c)
  : capacity(capacity),
    tree(c)
{
}



template <class T>
bool FixedSizePriorityQueue<T>::isEmpty() const
{
  return tree.isEmpty();
}



template <class T>
int FixedSizePriorityQueue<T>::getCapacity() const
{
  return capacity;
}



template <class T>
int FixedSizePriorityQueue<T>::getNumElements() const
{
  return tree.getCardinality();
}



template <class T>
bool FixedSizePriorityQueue<T>::insert(T t)
{
  if(tree.getCardinality()<capacity)
    {
      tree.insert(t);
      return true;
    }

  // Otherwise, the queue is at capacity

  // See if the new node is good enough to displace another
  T minElem=tree.minimum();
  if(tree.getComparator().less(t,minElem)) return false;
  tree.remove(minElem);
  tree.insert(t);
  return true;
}



template <class T>
T &FixedSizePriorityQueue<T>::peekMin() const
{
  return tree.minimum();
}



template <class T>
T &FixedSizePriorityQueue<T>::peekMax() const
{
  return tree.maximum();
}



template <class T>
FixedSizePriorityQueue<T>::iterator FixedSizePriorityQueue<T>::begin()
{
  return tree.begin();
}



template <class T>
FixedSizePriorityQueue<T>::iterator FixedSizePriorityQueue<T>::end()
{
  return tree.end();
}


#endif
