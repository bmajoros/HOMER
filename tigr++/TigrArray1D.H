/**************************************************************
TigrArray1D.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrArray1D_H
#define INCL_TigrArray1D_H
using namespace std;
#include <iostream>
#include <sstream>
#include "TigrStacktrace.H"
#include "TigrFile.H"

//#define DEBUG

template<class T>
class TigrArray1D
{
public:
  TigrArray1D(int size);
  TigrArray1D(const TigrArray1D &);
  virtual ~TigrArray1D();
  T &operator[](int);
  const T &operator[](int) const;
  int size() const;
  void printOn(ostream &os);
  void resize(int newSize);
  void setAllTo(T);
  void saveBytes(TigrFile &); // use only for basic data types!
  void loadBytes(TigrFile &); // use only for basic data types!
  TigrArray1D<T> &operator=(const TigrArray1D<T> &);
private:
  int n;
  T *array;
  void checkIndex(int index) const;
};

typedef TigrArray1D<double>   TigrDblArray1D;
typedef TigrArray1D<float>    TigrFloatArray1D;
typedef TigrArray1D<int>      TigrIntArray1D;
typedef TigrArray1D<long>     TigrLongArray1D;
typedef TigrArray1D<unsigned> TigrUnsignedArray1D;
typedef TigrArray1D<short>    TigrShortArray1D;
typedef TigrArray1D<char>     TigrCharArray1D;

template<class T>
ostream &operator<<(ostream &os,TigrArray1D<T> &a);



template<class T>
TigrArray1D<T>::TigrArray1D(const TigrArray1D &other)
  : array(new T[other.n]), n(other.n)
{
  for(int i=0 ; i<n ; ++i)
    array[i]=other.array[i];
}



template<class T>
TigrArray1D<T>::TigrArray1D(int ni)
  : n(ni), array(ni?new T[ni]:NULL)
{
  // ctor
}



template<class T>
TigrArray1D<T>::~TigrArray1D()
{
  delete [] array;
}



template<class T>
ostream &operator<<(ostream &os,TigrArray1D<T> &a)
{
  a.printOn(os);
  return os;
}



template<class T>
T &TigrArray1D<T>::operator[](int i)
{
#ifdef DEBUG
  checkIndex(i);
#endif
  return array[i];
}



template<class T>
const T &TigrArray1D<T>::operator[](int i) const
{
#ifdef DEBUG
  checkIndex(i);
#endif
  return array[i];
}



template<class T>
int TigrArray1D<T>::size() const
{
  return n;
}



template<class T>
void TigrArray1D<T>::checkIndex(int index) const
{
  if(index<0 || index>=n)
    {
      stringstream os;
      os << "Invalid index (" << index << ") in TigrArray1D"
	 << ends;
      throw TigrStacktrace(os.str());
    }
}



template<class T>
void TigrArray1D<T>::printOn(ostream &os)
{
  int i;
  for(i=0 ; i<n ; ++i)
    os << array[i] << '\t';
}



template<class T>
void TigrArray1D<T>::resize(int newSize)
{
  delete [] array;
  array=(newSize ? new T[newSize] : NULL);
  n=newSize;
}



template<class T>
void TigrArray1D<T>::setAllTo(T t)
{
  int i;
  for(i=0 ; i<n ; ++i)
    array[i]=t;
}



template<class T>
void TigrArray1D<T>::saveBytes(TigrFile &file)
{
  file << n;
  int size=n*sizeof(T);
  file.write(size,static_cast<void*>(array));
}



template<class T>
void TigrArray1D<T>::loadBytes(TigrFile &file)
{
  file >> n;
  int size=n*sizeof(T);
  resize(n);
  file.read(size,static_cast<void*>(array));
}



template<class T>
TigrArray1D<T> &TigrArray1D<T>::operator=(const TigrArray1D<T> &other)
{
  if(n!=other.n)
    {
      delete [] array;
      n=other.n;
      array=new T[n];
    }
  for(int i=0 ; i<n ; ++i) 
    array[i]=other.array[i];
  return *this;
}






#endif
