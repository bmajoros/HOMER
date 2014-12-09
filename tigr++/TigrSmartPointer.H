/****************************************************************
 TigrSmartPointer.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_TigrSmartPointer_H
#define INCL_TigrSmartPointer_H
#include <iostream>
using namespace std;

/****************************************************************
                   HOW TO USE SMART POINTERS

 Very simple: 

   (1) The class that will be pointed to must either be
       derived from class PointedTo (below) or have a 
       "referenceCount" integer attribute.  

   (2) Be sure to #define USE_SMART_POINTERS before including this
       header file.  (Then you can easily disable smart-pointing 
       simply by commenting out that #define)

   (3) Declare your pointers as:
            TigrSmartPointer<SomeClassName> ptr;
       rather than:
            SomeClassName *ptr;

*****************************************************************/




/****************************************************************
                         base class PointedTo
 */
struct PointedTo
{
#ifdef USE_SMART_POINTERS
  int referenceCount;
  PointedTo() : referenceCount(0) {}
#endif
};





/****************************************************************
                       class TigrSmartPointer
 */
template<class T>
class TigrSmartPointer
{
  T *p;
public:
  //-----------------------------------------------------
  TigrSmartPointer(T *p=NULL) 
    : p(p) 
  {
#ifdef USE_SMART_POINTERS
    if(p) ++p->referenceCount;
#endif
  }
  //-----------------------------------------------------
  TigrSmartPointer(const TigrSmartPointer<T> &o) 
    : p(o.p) 
  {
#ifdef USE_SMART_POINTERS
    if(p) ++p->referenceCount;
#endif
  }
  //-----------------------------------------------------
#ifdef USE_SMART_POINTERS
  ~TigrSmartPointer() 
  {
    if(p)
      if(--p->referenceCount==0) delete p;
  }
#endif
  //-----------------------------------------------------
  TigrSmartPointer<T> &operator=(const TigrSmartPointer<T> &o)
  {
#ifdef USE_SMART_POINTERS
    if(p==o.p) return *this;
    if(p) 
      if(--p->referenceCount==0) delete p;
#endif
    p=o.p;
#ifdef USE_SMART_POINTERS
    if(p) ++p->referenceCount; 
#endif
    return *this;
  }
  //-----------------------------------------------------
  TigrSmartPointer<T> &operator=(T *t)
  {
#ifdef USE_SMART_POINTERS
    if(p==t) return *this;
    if(p) 
      if(--p->referenceCount==0) delete p;
#endif
    p=t;
#ifdef USE_SMART_POINTERS
    if(p) ++p->referenceCount; 
#endif
    return *this;
  }
  //-----------------------------------------------------
  T *operator ->() {return p;}
  //-----------------------------------------------------
  const T *operator ->() const {return p;}
  //-----------------------------------------------------
#ifndef USE_SMART_POINTERS
  T &operator*() {return *p;}
#endif
  //-----------------------------------------------------
#ifndef USE_SMART_POINTERS
  const T &operator*() const {return *p;}
#endif
  //-----------------------------------------------------
#ifndef USE_SMART_POINTERS
  operator T*() const {return p;}
#endif
  //-----------------------------------------------------
#ifndef USE_SMART_POINTERS
  //operator const T*() const {return p;}
#endif
  //-----------------------------------------------------
  bool operator<(const TigrSmartPointer<T> &t) const {return p<t.p;}
  //-----------------------------------------------------
  bool operator==(const TigrSmartPointer<T> &t) const {return p==t.p;}
  //-----------------------------------------------------
  bool operator!=(const TigrSmartPointer<T> &t) const {return p!=t.p;}
  //-----------------------------------------------------
  bool operator>(const TigrSmartPointer<T> &t) const {return p>t.p;}
  //-----------------------------------------------------
  operator bool() {return p!=NULL;}
  //-----------------------------------------------------
  void printOn(ostream &os) const {os<<p;}
};

template<class T>
ostream &operator<<(ostream &os,const TigrSmartPointer<T> &p)
{
  p.printOn(os);
  return os;
}

#endif
