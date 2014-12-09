/****************************************************************
 Environment.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_Environment_H
#define INCL_Environment_H
using namespace std;
#include "TigrString.H"

class Environment
{
public:
  static TigrString lookup(const TigrString &);
};


#endif
