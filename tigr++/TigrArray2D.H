/**************************************************************
TigrArray2D.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrArray2D_H
#define INCL_TigrArray2D_H
using namespace std;
#include <iostream>
#include "TigrString.H"
#include "TigrStacktrace.H"
#include "TigrFile.H"

template<class T>
class TigrArray2D
{
public:
template<class X>
  class RowIn2DArray
  {
  public:
    RowIn2DArray(TigrArray2D<X> &,int);
    X &operator[](int);
    const X &operator[](int) const;
  private:
    int firstInd;
    TigrArray2D<X> &theArray;
  };
  TigrArray2D(int,int);
  TigrArray2D(const TigrArray2D &);
  virtual ~TigrArray2D();

RowIn2DArray<T> operator[](int i)
{
  return RowIn2DArray<T>(*this,i);
}

const RowIn2DArray<T> operator[](int i) const
{
  return RowIn2DArray<T>(const_cast<TigrArray2D<T>&>(*this),i);
}
//RowIn2DArray<T> operator[](int);
//const RowIn2DArray<T> operator[](int) const;
  T &index(int,int);
  int getFirstDim() const;
  int getSecondDim() const;
  void printOn(ostream &) const;
  void readFrom(istream &);
  void setAllTo(T);
  void resize(int,int);
  void saveBytes(TigrFile &); // use only for basic data types!
  void loadBytes(TigrFile &); // use only for basic data types!
private:
  int dim1, dim2;
  T *theArray;
};

template<class T>
ostream &operator<<(ostream &,const TigrArray2D<T> &);

template<class T>
istream &operator>>(istream &,TigrArray2D<T> &);


typedef TigrArray2D<double>   TigrDblArray2D;
typedef TigrArray2D<float>    TigrFloatArray2D;
typedef TigrArray2D<int>      TigrIntArray2D;
typedef TigrArray2D<long>     TigrLongArray2D;
typedef TigrArray2D<unsigned> TigrUnsignedArray2D;
typedef TigrArray2D<short>    TigrShortArray2D;
typedef TigrArray2D<char>     TigrCharArray2D;



template<class T>
TigrArray2D<T>::TigrArray2D(const TigrArray2D &other)
  : dim1(other.dim1), dim2(other.dim2), 
    theArray(new T[other.dim1*other.dim2])
{
  //### THIS IS NOT IMPLEMENTED OPTIMALLY -- BILL, FIX THIS...

  T *dest=theArray, *src=other.theArray;
  for(int i=0 ; i<dim1 ; ++i)
    for(int j=0 ; j<dim2 ; ++j)
      theArray[i*dim2+j]=other.theArray[i*dim2+j];
}



template<class T>
TigrArray2D<T>::TigrArray2D(int d1,int d2)
  : dim1(d1), dim2(d2), theArray(new T[d1*d2])
{
  // ctor
}



template<class T>
TigrArray2D<T>::~TigrArray2D()
{
  // dtor
  delete [] theArray;
}



template<class T>
ostream &operator<<(ostream &os,const TigrArray2D<T> &theArray)
{
  theArray.printOn(os);
  return os;
}



template<class T>
template<class X>
X &TigrArray2D<T>::RowIn2DArray<X>::operator[](int secondIndex)
{
  return theArray.index(firstInd,secondIndex);
}



template<class T>
template<class X>
const X &TigrArray2D<T>::RowIn2DArray<X>::operator[](int secondIndex) const
{
  return theArray.index(firstInd,secondIndex);
}



/*template<class T>
template<class X>
TigrArray2D<T>::RowIn2DArray<X> TigrArray2D<T>::operator[](int i)
{
  return RowIn2DArray<T>(*this,i);
}
*/

 /*
template<class T>
template<class X>
const TigrArray2D<T>::RowIn2DArray<X> TigrArray2D<T>::operator[](int i) 
  const
{
  return RowIn2DArray<T>(const_cast<TigrArray2D<T>&>(*this),i);
}
 */


template<class T>
template<class X>
TigrArray2D<T>::RowIn2DArray<X>::RowIn2DArray(TigrArray2D<X> &theArray,
					      int i)
  : theArray(theArray), firstInd(i)
{
  // ctor
}



template<class T>
T &TigrArray2D<T>::index(int i,int j)
{
#ifdef DEBUG
  if(i<0 || i>=dim1 || j<0 || j>=dim2)
    throw TigrStacktrace(TigrString("index out of bounds in TigrArray2D<>::index(")+i+","+j+") dim1="+dim1+" dim2="+dim2+")");
#endif

  int ind=i*dim2+j;
  return theArray[ind];
}



template<class T>
int TigrArray2D<T>::getFirstDim() const
{
  return dim1;
}



template<class T>
int TigrArray2D<T>::getSecondDim() const
{
  return dim2;
}



template<class T>
void TigrArray2D<T>::setAllTo(T f)
{
  int i, j;
  for(i=0 ; i<dim1 ; ++i)
    for(j=0 ; j<dim2 ; ++j)
      index(i,j)=f;
}


template<class T>
istream &operator>>(istream &is,TigrArray2D<T> &array)
{
  array.readFrom(is);
  return is;
}


template<class T>
void TigrArray2D<T>::printOn(ostream &os) const
{
  int i, j;
  for(j=0 ; j<dim2 ; ++j)
    {
      for(i=0 ; i<dim1 ; ++i)
	os << const_cast<TigrArray2D*>(this)->index(i,j) << '\t';
      os << endl;
    }
}



template<class T>
void TigrArray2D<T>::readFrom(istream &is)
{
  TigrString buf;
  for(int j=0 ; j<dim2 ; ++j)
    for(int i=0 ; i<dim1 ; ++i)

      {
	is >> buf;
	index(i,j)=(T) buf;
      }  
}


template<class T>
void TigrArray2D<T>::resize(int x,int y)
{
  delete [] theArray;
  dim1=x;
  dim2=y;
  theArray=new T[dim1*dim2];
}



template<class T>
void TigrArray2D<T>::saveBytes(TigrFile &file)
{
  file << dim1 << dim2;
  int size=sizeof(T)*dim1*dim2;
  file.write(size,static_cast<void*>(theArray));
}



template<class T>
void TigrArray2D<T>::loadBytes(TigrFile &file)
{
  file >> dim1 >> dim2;
  resize(dim1,dim2);
  int size=sizeof(T)*dim1*dim2;
  file.read(size,static_cast<void*>(theArray));
}


#endif
