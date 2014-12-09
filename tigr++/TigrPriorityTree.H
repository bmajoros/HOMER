/**************************************************************
TigrPriorityTree.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrPriorityTree_H
#define INCL_TigrPriorityTree_H
#include "TigrFastBinTree.H"
using namespace std;

class TigrPriorityTree
{
public:
  TigrPriorityTree();
  TigrFastBinNode *extractMin();
  TigrFastBinNode *peekMin();
  void insert(TigrFastBinNode *);
private:
  TigrFastBinTree tree;
};

#endif
