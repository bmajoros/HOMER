/****************************************************************
 TigrIterator.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/

#ifndef INCL_TigrIterator_H
#define INCL_TigrIterator_H
using namespace std;

template <class T>
class TigrIterator
{
public:
  virtual void operator++()=0;
  virtual void operator++(int)=0;
  virtual T &operator*()=0;
  virtual bool operator==(const TigrIterator<T> &) const=0;
  virtual bool operator!=(const TigrIterator<T> &i) const
    {return !(i==*this);}
  virtual TigrIterator<T> &clone()=0;
};


#endif
