/**************************************************************
TigrHashpjw.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrHashpjw_H
#define INCL_TigrHashpjw_H
#include "TigrString.H"
using namespace std;

unsigned TigrHashpjw(TigrString,unsigned tableSize);

unsigned TigrHashpjw(const char *,unsigned tableSize);

unsigned TigrHashpjw(const TigrString &,unsigned tableSize);

#endif
