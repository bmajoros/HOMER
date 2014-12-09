/**************************************************************
TigrVectorSorter.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrVectorSorter_H
#define INCL_TigrVectorSorter_H
#include <math.h>
#include "TigrVector.H"
#include "TigrComparator.H"
#include "TigrVectorPrinter.H"
#include "TigrRandom.H"
using namespace std;

/****************************************************************
  TigrVectorSorter -- STANDARD QUICK-SORT ALGORITHM

  
*****************************************************************/
template <class T>
class TigrVectorSorter
{
public:
  TigrVectorSorter(TigrVector<T> &,TigrComparator<T> &);
  TigrVector<T> *sortAscending();
  TigrVector<T> *sortDescending();  
  TigrVector<int> *sortAscendByIndex();
  TigrVector<int> *sortDescendByIndex();
  void sortAscendInPlace();
  void sortDescendInPlace();
  void applyIndexMap(TigrVector<int> &indices,TigrVector<T> &from,
		     TigrVector<T> &to);
  void applyIndexMap(TigrVector<int> &indices,TigrVector<T> &fromAndTo);
private:
  int numElements;
  TigrVector<T> &theVector;
  TigrComparator<T> &comp;

  TigrVector<T> *quickSort();
  TigrVector<int> *sortByIndex();
  bool (TigrComparator<T>::*ptrToMmbrGreater)(T &,T &);
  bool (TigrComparator<T>::*ptrToMmbrLess)(T &,T &);
  int partition(TigrVector<T> &,int begin,int end);
  int partitionIndices(TigrVector<int> &indices,TigrVector<T> &,
		       int begin,int end);
  void initIndexArray(TigrVector<int> &);
  void quickSort(TigrVector<T> &,int begin,int end);
  void quickSortIndices(TigrVector<int> &indices,TigrVector<T> &,
			int begin,int end);
  void sortInPlace();
  void swap(TigrVector<T> &,int,int);
  void swapIndices(TigrVector<int> &,int,int);
};



template <class T>
TigrVectorSorter<T>::TigrVectorSorter(TigrVector<T> &v,TigrComparator<T> &c)
  : theVector(v), numElements(v.size()), comp(c)
{
}



template <class T>
TigrVector<T> *TigrVectorSorter<T>::quickSort()
{
  TigrVector<T> &sorted=*new TigrVector<T>(numElements);
  for(int i=0 ; i<numElements ; ++i) sorted[i]=theVector[i];
  quickSort(sorted,0,numElements-1);
  return &sorted;
}



template <class T>
TigrVector<T> *TigrVectorSorter<T>::sortAscending()
{
  ptrToMmbrLess=&TigrComparator<T>::less;
  ptrToMmbrGreater=&TigrComparator<T>::greater;
  return quickSort();
}



template <class T>
TigrVector<T> *TigrVectorSorter<T>::sortDescending()
{
  ptrToMmbrLess=&TigrComparator<T>::greater;
  ptrToMmbrGreater=&TigrComparator<T>::less;
  return quickSort();
}



template <class T>
TigrVector<int> *TigrVectorSorter<T>::sortByIndex()
{
  TigrVector<int> &indices=*new TigrVector<int>;
  initIndexArray(indices);
  quickSortIndices(indices,theVector,0,numElements-1);
  return &indices;
}



template <class T>
TigrVector<int> *TigrVectorSorter<T>::sortAscendByIndex()
{
  ptrToMmbrLess=&TigrComparator<T>::less;
  ptrToMmbrGreater=&TigrComparator<T>::greater;
  return sortByIndex();
}



template <class T>
TigrVector<int> *TigrVectorSorter<T>::sortDescendByIndex()
{
  ptrToMmbrLess=&TigrComparator<T>::greater;
  ptrToMmbrGreater=&TigrComparator<T>::less;
  return sortByIndex();
}



template <class T>
int TigrVectorSorter<T>::partition(TigrVector<T> &sorted,
				   int begin,int end)
{
   int pivotIndex=RandomNumber(end-begin+1)+begin;
   T temp=sorted[pivotIndex];
   swap(sorted,pivotIndex,begin);
   int i=begin-1;
   int j=end+1;
   while(true) 
     {
       do --j;
       while((comp.*ptrToMmbrGreater)(sorted[j],temp));
       do ++i;
       while((comp.*ptrToMmbrLess)(sorted[i],temp));
       if(i>=j) break;
       else swap(sorted,i,j);
     }
   return j;
}



template <class T>
void TigrVectorSorter<T>::initIndexArray(TigrVector<int> &indices)
{
  for(int i=0 ; i<numElements ; ++i) indices.push_back(i);
}



template <class T>
void TigrVectorSorter<T>::quickSort(TigrVector<T> &sorted,
			int from,int to)
{
  if(from>=to) return;
  int pivot=partition(sorted,from,to);
  quickSort(sorted,from,pivot);
  quickSort(sorted,pivot+1,to);
}



template <class T>
void TigrVectorSorter<T>::sortInPlace()
{
  quickSort(theVector,0,numElements-1);
}



template <class T>
void TigrVectorSorter<T>::sortAscendInPlace()
{
  ptrToMmbrGreater=&TigrComparator<T>::greater;
  ptrToMmbrLess=&TigrComparator<T>::less;
  sortInPlace();
}



template <class T>
void TigrVectorSorter<T>::sortDescendInPlace()
{
  ptrToMmbrGreater=&TigrComparator<T>::less;
  ptrToMmbrLess=&TigrComparator<T>::greater;
  return sortInPlace();
}



template <class T>
void TigrVectorSorter<T>::swap(TigrVector<T> &data,int index1,int index2)
{
  T temp=data[index1];
  data[index1]=data[index2];
  data[index2]=temp;
}



template <class T>
void TigrVectorSorter<T>::swapIndices(TigrVector<int> &indices,
				      int index1,int index2)
{
  int temp=indices[index1];
  indices[index1]=indices[index2];
  indices[index2]=temp;
}



template <class T>
void TigrVectorSorter<T>::applyIndexMap(TigrVector<int> &indices,
					TigrVector<T> &from,
					TigrVector<T> &to)
{
  int n=from.size();
  for(int i=0 ; i<n ; ++i)
    to[i]=from[indices[i]];
}



template <class T>
void TigrVectorSorter<T>::applyIndexMap(TigrVector<int> &indices,
					TigrVector<T> &fromAndTo)
{
  int n=fromAndTo.size();
  TigrVector<T> temp(n);
  applyIndexMap(indices,fromAndTo,temp);
  for(int i=0 ; i<n ; ++i)
    fromAndTo[i]=temp[i];
}



template <class T>
void TigrVectorSorter<T>::quickSortIndices(TigrVector<int> &indices,
					   TigrVector<T> &data,
					   int from,int to)
{
  if(from>=to) return;
  int pivot=partitionIndices(indices,data,from,to);
  quickSortIndices(indices,data,from,pivot);
  quickSortIndices(indices,data,pivot+1,to);
}



template <class T>
int TigrVectorSorter<T>::partitionIndices(TigrVector<int> &indices,
					  TigrVector<T> &data,
					  int begin,int end)
{
   int pivotIndex=RandomNumber(end-begin+1)+begin;
   T temp=data[indices[pivotIndex]];
   swapIndices(indices,pivotIndex,begin);
   int i=begin-1;
   int j=end+1;
   while(true) 
     {
       do --j;
       while((comp.*ptrToMmbrGreater)(data[indices[j]],temp));
       do ++i;
       while((comp.*ptrToMmbrLess)(data[indices[i]],temp));
       if(i>=j) break;
       else swapIndices(indices,i,j);
     }
   return j;
}


#endif
