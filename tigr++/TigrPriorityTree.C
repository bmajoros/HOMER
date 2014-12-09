#include <iostream>
#include "TigrPriorityTree.H"
using namespace std;



TigrPriorityTree::TigrPriorityTree()
{
   
}



TigrFastBinNode *TigrPriorityTree::extractMin()
{
  TigrFastBinNode *m=dynamic_cast<TigrFastBinNode*>(tree.minimum());
  if(m) tree.deleteNode(m);
  return m;
}



TigrFastBinNode *TigrPriorityTree::peekMin()
{
  return dynamic_cast<TigrFastBinNode*>(tree.minimum());
}



void TigrPriorityTree::insert(TigrFastBinNode *node)
{
  tree.insert(node);
}
