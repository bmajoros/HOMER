/****************************************************************
 MemoryProfiler.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_MemoryProfiler_H
#define INCL_MemoryProfiler_H
using namespace std;
#include "TigrString.H"
#include <iostream>

class MemoryProfiler
{
public:
  static void report(const TigrString &label,ostream & =cout);
};


#endif
