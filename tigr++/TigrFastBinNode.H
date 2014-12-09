/**************************************************************
TigrFastBinNode.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrFastBinNode_H
#define INCL_TigrFastBinNode_H
#include "TigrBaseBinTreeNode.H"
using namespace std;

class TigrFastBinNode : public TigrBaseBinTreeNode
{
public:
  TigrFastBinNode *parent;
};

#endif
