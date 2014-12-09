/****************************************************************
 swap.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/

#ifndef INCL_swap_H
#define INCL_swap_H


template<class T>
inline void swap(T &a,T &b)
{
  T c=a;
  a=b;
  b=c;
}


#endif
