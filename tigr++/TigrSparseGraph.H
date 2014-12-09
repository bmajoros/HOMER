/**************************************************************
TigrSparseGraph.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrSparseGraph_H
#define INCL_TigrSparseGraph_H
#include "TigrVector.H"
#include "TigrSet.H"
using namespace std;

typedef unsigned VertexId;
typedef TigrSet<VertexId> Neighborhood;
class TigrSparseGraph
{
public:
  TigrSparseGraph(int numVertices=0);
  Neighborhood &getNeighborsOf(VertexId);
  VertexId addVertex();
  bool areAdjecent(VertexId,VertexId);
  int getDegree(VertexId);
  int getNumEdges();
  int getNumVertices();
  void addEdge(VertexId from,VertexId to);
  void removeEdge(VertexId from,VertexId to);
private:
  int numEdges;
  TigrVector<Neighborhood> neighborhoods;
};

#endif
