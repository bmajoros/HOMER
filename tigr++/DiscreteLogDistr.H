/****************************************************************
 DiscreteLogDistr.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_DiscreteLogDistr_H
#define INCL_DiscreteLogDistr_H
#include "DiscreteDistr.H"
#include "TigrArray1D.H"


template<class FLOAT_OR_DOUBLE>
class DiscreteLogDistr : public DiscreteDistr<FLOAT_OR_DOUBLE>
{
  virtual void initialize(const ContinuousDistr &);
public:
  DiscreteLogDistr(int);
  DiscreteLogDistr(int,const ContinuousDistr &);
};



template<class T>
DiscreteLogDistr<T>::DiscreteLogDistr(int maxValue)
  : DiscreteDistr<T>(maxValue)
{
}



template<class T>
DiscreteLogDistr<T>::DiscreteLogDistr(int maxValue,
			      const ContinuousDistr &d)
  : DiscreteDistr<T>(maxValue)
{
  initialize(d);
}



template<class T>
void DiscreteLogDistr<T>::initialize(const ContinuousDistr &d)
{
  mean=0;
  T total=0;
  int size=array.size();
  for(int i=0 ; i<size ; ++i)
    {
      T p=d.density(double(i));
      array[i]=p;
      total+=p;
    }
  for(int i=0 ; i<size ; ++i)
    {
      array[i]/=total;
      mean+=array[i]*i;
    }
  for(int i=0 ; i<size ; ++i)
    array[i]=log(array[i]);
  Pmean=d.density(double(mean));
}




#endif
