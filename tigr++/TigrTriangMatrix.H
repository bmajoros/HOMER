/**************************************************************
TigrTriangMatrix.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrTriangMatrix_H
#define INCL_TigrTriangMatrix_H
#include "TigrVector.H"
#include "TigrExceptions.H"
using namespace std;

template<class T>
class TigrTriangMatrix
{
public:
  TigrTriangMatrix(int n);
  T &operator()(int i,int j); 
  int getN() const;
private:
  int arraySize;
  int n;
  TigrVector<T> array;
  int computeIndex(int i,int j);
};



template<class T>
T &TigrTriangMatrix<T>::operator()(int i,int j)
{
  if(i>j)
    throw RootException("TigrTriangMatrix::operator(): i>j");
  int index=computeIndex(i,j);
  return array[index];
}



template<class T>
TigrTriangMatrix<T>::TigrTriangMatrix(int N)
  : array(((N+1)*N)/2), n(N)
{
  arraySize=array.size();
}



template<class T>
int TigrTriangMatrix<T>::computeIndex(int i,int j)
{
  return arraySize-(n-i-1)*(n-i)/2-n+j;
}



template<class T>
int TigrTriangMatrix<T>::getN() const
{
  return n;
}


#endif
