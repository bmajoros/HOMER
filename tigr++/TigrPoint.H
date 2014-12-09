/****************************************************************
 TigrPoint.H

 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/

#ifndef INCL_TigrPoint_H
#define INCL_TigrPoint_H
using namespace std;

template<class T>
struct TigrPoint
{
  T x, y;
  TigrPoint(T,T);
  TigrPoint() {}
};


template<class T>
TigrPoint<T>::TigrPoint(T x,T y)
{
  this->x=x;
  this->y=y;
}


#endif
