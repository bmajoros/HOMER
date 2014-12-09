/**************************************************************
TigrArray1DSorter.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrArray1DSorter_H
#define INCL_TigrArray1DSorter_H
#include <math.h>
#include "TigrArray1D.H"
#include "TigrComparator.H"
#include "TigrRandom.H"
using namespace std;

/****************************************************************
  TigrArray1DSorter -- STANDARD QUICK-SORT ALGORITHM

  
*****************************************************************/
template <class T>
class TigrArray1DSorter
{
public:
  TigrArray1DSorter(TigrArray1D<T> &,TigrComparator<T> &);
  TigrArray1D<T> *sortAscending();
  TigrArray1D<T> *sortDescending();  
  TigrArray1D<int> *sortAscendByIndex();
  TigrArray1D<int> *sortDescendByIndex();
  void sortAscendInPlace();
  void sortDescendInPlace();
private:
  int numElements;
  TigrArray1D<T> &theVector;
  TigrComparator<T> &comp;

  TigrArray1D<T> *quickSort();
  TigrArray1D<int> *sortByIndex();
  bool (TigrComparator<T>::*ptrToMmbrGreater)(T &,T &);
  bool (TigrComparator<T>::*ptrToMmbrLess)(T &,T &);
  int partition(TigrArray1D<T> &,int begin,int end);
  int partition(TigrArray1D<int> &,int begin,int end);
  void initIndexArray(TigrArray1D<int> &);
  void quickSort(TigrArray1D<T> &,int begin,int end);
  void quickSort(TigrArray1D<int> &,int begin,int end);
  void sortInPlace();
  void swap(TigrArray1D<T> &,int,int);
  void swap(TigrArray1D<int> &,int,int);
};



template <class T>
TigrArray1DSorter<T>::TigrArray1DSorter(TigrArray1D<T> &v,TigrComparator<T> &c)
  : theVector(v), numElements(v.size()), comp(c)
{
}



template <class T>
TigrArray1D<T> *TigrArray1DSorter<T>::quickSort()
{
  TigrArray1D<T> &sorted=*new TigrArray1D<T>(numElements);
  for(int i=0 ; i<numElements ; ++i) sorted[i]=theVector[i];
  quickSort(sorted,0,numElements-1);
  return &sorted;
}



template <class T>
TigrArray1D<T> *TigrArray1DSorter<T>::sortAscending()
{
  ptrToMmbrLess=&TigrComparator<T>::less;
  ptrToMmbrGreater=&TigrComparator<T>::greater;
  return quickSort();
}



template <class T>
TigrArray1D<T> *TigrArray1DSorter<T>::sortDescending()
{
  ptrToMmbrLess=&TigrComparator<T>::greater;
  ptrToMmbrGreater=&TigrComparator<T>::less;
  return quickSort();
}



template <class T>
TigrArray1D<int> *TigrArray1DSorter<T>::sortByIndex()
{
  TigrArray1D<int> &indices=*new TigrArray1D<int>;
  initIndexArray(indices);
  quickSort(indices,0,numElements-1);
  return &indices;
}



template <class T>
TigrArray1D<int> *TigrArray1DSorter<T>::sortAscendByIndex()
{
  ptrToMmbrLess=&TigrComparator<T>::less;
  ptrToMmbrGreater=&TigrComparator<T>::greater;
  return sortByIndex();
}



template <class T>
TigrArray1D<int> *TigrArray1DSorter<T>::sortDescendByIndex()
{
  ptrToMmbrLess=&TigrComparator<T>::greater;
  ptrToMmbrGreater=&TigrComparator<T>::less;
  return sortByIndex();
}



template <class T>
int TigrArray1DSorter<T>::partition(TigrArray1D<T> &sorted,
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
int TigrArray1DSorter<T>::partition(TigrArray1D<int> &indices,
				   int begin,int end)
{
   int pivotIndex=RandomNumber(end-begin+1)+begin;
   int tempIndex=indices[pivotIndex];
   T &temp=theVector[tempIndex];
   swap(indices,pivotIndex,begin);
   int i=begin-1;
   int j=end+1;
   while(true) 
     {
       do --j;
       while((comp.*ptrToMmbrGreater)(theVector[indices[j]],temp));
       do ++i;
       while((comp.*ptrToMmbrLess)(theVector[indices[i]],temp));
       if(i>=j) break;
       else swap(indices,i,j);
     }
   return j;
}



template <class T>
void TigrArray1DSorter<T>::initIndexArray(TigrArray1D<int> &indices)
{
  for(int i=0 ; i<numElements ; ++i) indices.push_back(i);
}



template <class T>
void TigrArray1DSorter<T>::quickSort(TigrArray1D<T> &sorted,
			int from,int to)
{
  if(from>=to) return;
  int pivot=partition(sorted,from,to);
  quickSort(sorted,from,pivot);
  quickSort(sorted,pivot+1,to);
}



template <class T>
void TigrArray1DSorter<T>::quickSort(TigrArray1D<int> &indices,int from,int to)
{
  if(from>=to) return;
  int pivot=partition(indices,from,to);
  quickSort(indices,from,pivot);
  quickSort(indices,pivot+1,to);
}



template <class T>
void TigrArray1DSorter<T>::sortInPlace()
{
  quickSort(theVector,0,numElements-1);
}



template <class T>
void TigrArray1DSorter<T>::sortAscendInPlace()
{
  ptrToMmbrGreater=&TigrComparator<T>::greater;
  ptrToMmbrLess=&TigrComparator<T>::less;
  sortInPlace();
}



template <class T>
void TigrArray1DSorter<T>::sortDescendInPlace()
{
  ptrToMmbrGreater=&TigrComparator<T>::less;
  ptrToMmbrLess=&TigrComparator<T>::greater;
  return sortInPlace();
}



template <class T>
void TigrArray1DSorter<T>::swap(TigrArray1D<T> &indices,int index1,int index2)
{
  T temp=indices[index1];
  indices[index1]=indices[index2];
  indices[index2]=temp;
}



template <class T>
void TigrArray1DSorter<T>::swap(TigrArray1D<int> &indices,int index1,int index2)
{
  int temp=indices[index1];
  indices[index1]=indices[index2];
  indices[index2]=temp;
}


#endif
