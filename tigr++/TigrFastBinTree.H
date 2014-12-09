/**************************************************************
TigrFastBinTree.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrFastBinTree_H
#define INCL_TigrFastBinTree_H
#include "TigrBinaryTree.H"
using namespace std;

/*============================================================
  TigrFastBinTree -- a binary tree template in which each node
  has a pointer to its left and right children, and ALSO a
  pointer to its parent.  Thus, finding a node's parent does NOT
  require searching down from the root, and is therefore very
  fast.  However, more memory is used by that extra pointer,
  and extra time for each insertion to initialize the parent 
  pointer.  This class also allows duplicate elements to be stored
  in the tree (if desired).
  ============================================================*/

template<class T>
struct TigrFastBinNode : public TigrBaseBinTreeNode<T>
{
  TigrFastBinNode<T> *parent;
  TigrFastBinNode() : parent(NULL) {}
};


template <class T>
class TigrFastBinTree : public TigrBinaryTree<T>
{
public:
  TigrFastBinTree(TigrComparator<T>& = *new DirectComparator<T>());
  virtual TigrBaseBinTreeNode<T> *getParent(TigrBaseBinTreeNode<T> *);
  virtual bool insert(T);
  virtual bool insert(TigrBaseBinTreeNode<T> *);
  virtual bool naiveInsert(TigrFastBinNode<T> *);
  virtual void deleteNode(TigrBaseBinTreeNode<T> *);
  void enableDuplicates();
private:
  bool allowDuplicates;
};



template <class T>
TigrFastBinTree<T>::TigrFastBinTree(TigrComparator<T> &c) 
  : allowDuplicates(false),
    TigrBinaryTree<T>(c)
{
}



template <class T>
TigrBaseBinTreeNode<T> *TigrFastBinTree<T>::getParent(
   TigrBaseBinTreeNode<T> *thisNode)
{
  TigrFastBinNode<T> *fastNode=static_cast<TigrFastBinNode<T>*>(thisNode);
  return fastNode->parent;
}



template <class T>
bool TigrFastBinTree<T>::insert(TigrBaseBinTreeNode<T> *newNode) 
{ 
  return naiveInsert(static_cast<TigrFastBinNode<T>*>(newNode)); 
}



template <class T>
bool TigrFastBinTree<T>::naiveInsert(TigrFastBinNode<T> *newNode)
{
  if(!root)
    {	
      root=newNode;
      newNode->left=newNode->right=newNode->parent=NULL;
      cardinality=1;
      return true;
    }
   
  TigrFastBinNode<T> *thisNode=static_cast<TigrFastBinNode<T>*>(root);
  while(true) 
    {
      // CASE -1 : LESS
      if(comparator.less(newNode->element,thisNode->element))
	{
	  if(!thisNode->left)
	    {
	      thisNode->left=newNode;
	      newNode->left=newNode->right=NULL;
	      newNode->parent=thisNode;
	      ++cardinality;
	      return true;
	    }
	  thisNode=static_cast<TigrFastBinNode<T>*>(thisNode->left);  
	  continue;
	}

      // CASE 0 : EQUAL
      else if(comparator.equal(newNode->element,thisNode->element))
	{
          if(!allowDuplicates)
	    return false;
	  // otherwise, fall through to CASE 1 below
	}

      // CASE 1 : GREATER
      if(!thisNode->right)
	{
	  thisNode->right=newNode;
	  newNode->left=newNode->right=NULL;
	  newNode->parent=thisNode;
	  ++cardinality;
	  return true;
	}
      
      thisNode=static_cast<TigrFastBinNode<T>*>(thisNode->right);
      continue;
    }
}



template <class T>
void TigrFastBinTree<T>::deleteNode(TigrBaseBinTreeNode<T> *baseNode)
{ 
  TigrFastBinNode<T> *thisNode=static_cast<TigrFastBinNode<T>*>(baseNode);
  int numChildren=0;
  TigrFastBinNode<T> *child, *mySuccessor, *hisParent;
  if(thisNode->left) numChildren++;
  if(thisNode->right) numChildren++;
   
  switch(numChildren) 
    {
    case 0:
      if(thisNode->parent)
	if(thisNode->parent->left==thisNode)
	  thisNode->parent->left=NULL;
	else
	  thisNode->parent->right=NULL;
      else  
	root=NULL;
      thisNode->parent=NULL;
      break;
      
    case 1:  
      if(thisNode->left)
	{
	  child=static_cast<TigrFastBinNode<T>*>(thisNode->left);
	  thisNode->left=NULL;
	}
      else
	{
	  child=static_cast<TigrFastBinNode<T>*>(thisNode->right);
	  thisNode->right=NULL;
	}
       
      if(thisNode->parent)
	{
	  if(thisNode->parent->left==thisNode)
	    thisNode->parent->left=child;
	  else
	    thisNode->parent->right=child;
	  child->parent=thisNode->parent;
	}
      else  
	{
	  root=child;
	  child->parent=NULL;
	}
      thisNode->parent=NULL;
      break;
      
    case 2:
      mySuccessor=static_cast<TigrFastBinNode<T>*>(thisNode->right);
      hisParent=NULL;
      while(mySuccessor->left)
	{
	  hisParent=mySuccessor;
	  mySuccessor=static_cast<TigrFastBinNode<T>*>(mySuccessor->left);
	}
      
      if(hisParent)  
	{
	  hisParent->left=mySuccessor->right;
	  if(mySuccessor->right)
	    static_cast<TigrFastBinNode<T>*>(mySuccessor->right)->parent=
	      hisParent;

	  mySuccessor->left=thisNode->left;
	  static_cast<TigrFastBinNode<T>*>(thisNode->left)->parent=mySuccessor;
	  
	  mySuccessor->right=thisNode->right;
	  static_cast<TigrFastBinNode<T>*>(thisNode->right)->parent=
	    mySuccessor;
	  
	  if(thisNode->parent)
	    {
	      if(thisNode->parent->left==thisNode)
		thisNode->parent->left=mySuccessor;
	      else
		thisNode->parent->right=mySuccessor;
	      mySuccessor->parent=thisNode->parent;
	    }
	  else  
	    {
	      root=mySuccessor;
	      mySuccessor->parent=NULL;
	    }
	  thisNode->left=thisNode->right=thisNode->parent=NULL;
	}
      
      else  
	{         	
	  mySuccessor->left=thisNode->left;
	  static_cast<TigrFastBinNode<T>*>(thisNode->left)->parent=
	    mySuccessor;
	  
	  if(thisNode->parent)
	    {
	      if(thisNode->parent->left==thisNode)
		thisNode->parent->left=mySuccessor;
	      else
		thisNode->parent->right=mySuccessor;
	      mySuccessor->parent=thisNode->parent;
	    }
	  else  
	    {
	      root=mySuccessor;
	      mySuccessor->parent=NULL;
	    }
	  
	  thisNode->left=thisNode->right=thisNode->parent=NULL;
	}
      break;
    }

  --cardinality;
  delete thisNode;
}



template <class T>
void TigrFastBinTree<T>::enableDuplicates() 
{ 
  allowDuplicates=true; 
}



template <class T>
bool TigrFastBinTree<T>::insert(T t)
{
  TigrFastBinNode<T> *n=new TigrFastBinNode<T>;
  n->element=t;
  if(insert(n)) return true;
  delete n;
  return false;
}


#endif
