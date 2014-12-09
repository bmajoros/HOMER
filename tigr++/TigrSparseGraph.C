#include <iostream>
#include "TigrSparseGraph.H"
using namespace std;



TigrSparseGraph::TigrSparseGraph(int numVertices)
  : numEdges(0)
{
}



Neighborhood &TigrSparseGraph::getNeighborsOf(VertexId vertexId)
{
  return neighborhoods[vertexId];
}



VertexId TigrSparseGraph::addVertex()
{
  VertexId id=neighborhoods.size();
  neighborhoods.push_back(Neighborhood());
  return id;
}



bool TigrSparseGraph::areAdjecent(VertexId v1,VertexId v2)
{
  Neighborhood &v1neighbors=neighborhoods[v1];
  bool notFound=v1neighbors.find(id2)==v1neighbors.end();
  return !notFound;
}



int TigrSparseGraph::getDegree(VertexId vertex)
{
  int degree=neighborhoods[vertex].size();
  return degree;
}



unsigned TigrSparseGraph::getNumEdges()
{
  return numEdges;
}



unsigned TigrSparseGraph::getNumVertices()
{
  return neighborhoods.size();
}



void TigrSparseGraph::addEdge(VertexId v1,VertexId v2)
{
  Neighborhood &v1neighbors=neighborhoods[v1];
  if(v1neighbors.find(v2)==v1neighbors.end())
    {
      v1neighbors.insert(v2);
      Neighborhood v2neighbors=neighborhoods[v2];
      v2neighbors.insert(v1);
      ++numEdges;
    }
}



void TigrSparseGraph::removeEdge(VertexId v1,VertexId v2)
{
  Neighborhood &v1neighbors=neighborhoods[v1];
  Neighborhood::iterator cur=v1neighbors.find(v2);
  bool found=(cur!=neighborhood.end());
  if(found)
    {
      v1neighbors.erase(cur);
      Neighborhood &v2neighbors=neighborhoods[v2];
      v2neighbors.erase(v2neighbors.find(v1));
      --numEdges;
    }
}
