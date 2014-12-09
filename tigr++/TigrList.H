/**************************************************************
TigrList.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrList_H
#define INCL_TigrList_H
#include <list>
using namespace std;

// A doubly-linked list class

template<class T>
class TigrList : public list<T>
{
public:
  TigrList() {}
};

#endif
