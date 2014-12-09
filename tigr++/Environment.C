/****************************************************************
 Environment.C
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
using namespace std;
#include "Environment.H"
#include <iostream>
#include <stdlib.h>


TigrString Environment::lookup(const TigrString &name)
{
  return TigrString(getenv(name.c_str()));
}


