/**************************************************************
TigrBinaryTree.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrBinaryTree_H
#define INCL_TigrBinaryTree_H
using namespace std;
#include "TigrComparator.H"

/*============================================================
  TigrBinaryTree -- a binary tree template in which each node
  has a pointer to its left and right children, but lacks a 
  pointer to its parent.  Thus, finding a node's parent requires
  searching down from the root.  If this is unacceptably slow,
  use the TigrFastBinaryTree instead (each node has a pointer to
  its parent).
  ============================================================*/

template<class T>
struct TigrBaseBinTreeNode
{
  TigrBaseBinTreeNode<T> *left, *right;
  T element;
  TigrBaseBinTreeNode() :left(NULL),right(NULL){}
};


template<class T>
class TigrBinaryTree
{
public:
  class TigrBinTreeIter
  {
    friend TigrBinaryTree<T>;
  public:
    enum IterationType
      {
	INORDER,   // tested
	POSTORDER, // tested
	PREORDER   // ### not tested
      };
    TigrBinTreeIter(TigrBinaryTree<T> *,IterationType=INORDER);
    T &operator*();
    operator bool();
    TigrBinTreeIter<T> &operator++();
    bool operator==(const TigrBinTreeIter<T> &) const;
  protected:
    enum FromType
      { 
	FROM_LEFT, 
	FROM_RIGHT, 
	FROM_PARENT, 
	FROM_NOWHERE 
      };
    TigrBinTreeIter(TigrBinaryTree<T> *,TigrBaseBinTreeNode<T> *thisNode,
		    TigrBaseBinTreeNode<T> *nextNode,FromType,
		    IterationType);
    virtual TigrBaseBinTreeNode<T> *inOrder();
    virtual TigrBaseBinTreeNode<T> *postOrder();
    virtual TigrBaseBinTreeNode<T> *preOrder();
    virtual void reset();  
    IterationType iterationType;
    FromType fromWhere;
    TigrBinaryTree<T> *myTree;
    TigrBaseBinTreeNode<T> *thisNode, *nextNode;
  };
  typedef TigrBinTreeIter<T> iterator;
  
  // ### The following constructor will cause a memory leak if you
  //     use a large number of these binary trees with the default 
  //     argument:
  TigrBinaryTree(TigrComparator<T>& = *new DirectComparator<T>());
  virtual ~TigrBinaryTree();                   // tested
  virtual bool insert(T);                      // tested
  virtual T &maximum() const;                  // tested
  virtual T &minimum() const;                  // tested
  virtual unsigned getCardinality();           
  bool isEmpty() const;                        // tested
  virtual void purge();                        // tested
  iterator begin();                            // tested
  iterator end();                              // tested
  void remove(T);                              // tested
  bool getPredecessor(T of,T &is);             // tested
  bool getSuccessor(T of,T &is);               // tested
  TigrComparator<T> &getComparator();
protected:
  friend TigrBinTreeIter;
  TigrComparator<T> &comparator;
  unsigned cardinality;
  TigrBaseBinTreeNode<T> *root;
  virtual TigrBaseBinTreeNode<T> *subTreeMax(TigrBaseBinTreeNode<T> *) const;
  virtual TigrBaseBinTreeNode<T> *subTreeMin(TigrBaseBinTreeNode<T> *) const;
  virtual TigrBaseBinTreeNode<T> *getParent(TigrBaseBinTreeNode<T> *);
  virtual TigrBaseBinTreeNode<T> *getRoot();
  virtual TigrBaseBinTreeNode<T> *predecessor(TigrBaseBinTreeNode<T> *);
  virtual TigrBaseBinTreeNode<T> *successor(TigrBaseBinTreeNode<T> *);
  virtual bool insert(TigrBaseBinTreeNode<T> *);
  virtual void deleteNode(TigrBaseBinTreeNode<T> *);
  TigrBaseBinTreeNode<T> *find(T);
};



template<class T>
TigrBinaryTree<T>::TigrBinaryTree(TigrComparator<T> &c) 
  : root(NULL),
    cardinality(0),
    comparator(c)
{
}



template<class T>
TigrBinaryTree<T>::~TigrBinaryTree()
{
  purge();
}



template<class T>
TigrComparator<T> &TigrBinaryTree<T>::getComparator()
{
  return comparator;
}



template<class T>
TigrBaseBinTreeNode<T> *TigrBinaryTree<T>::getParent(TigrBaseBinTreeNode<T> *
						  thisNode)
{
  if(thisNode==root) return NULL;

  TigrBaseBinTreeNode<T> *parent=root;
  while(true)
    {
      if(comparator.less(thisNode->element,parent->element))
	{
	  if(parent->left==thisNode) 
	    return parent;
	  parent=parent->left;
	  continue;
	}
      else
	{
	  if(parent->right==thisNode) 
	    return parent;
	  parent=parent->right;
	  continue;
	}
    }
}



template<class T>
TigrBaseBinTreeNode<T> *TigrBinaryTree<T>::getRoot()
{
  return root;
}



template<class T>
T &TigrBinaryTree<T>::maximum() const
{
  return subTreeMax(root)->element;
}



template<class T>
T &TigrBinaryTree<T>::minimum() const
{
  return subTreeMin(root)->element;
}



template<class T>
TigrBaseBinTreeNode<T> *TigrBinaryTree<T>::predecessor(TigrBaseBinTreeNode<T> *
						       ofNode)
{
  TigrBaseBinTreeNode<T> *thisNode;
   
  if(ofNode->left) 
    return subTreeMax(ofNode->left);
   
  thisNode=getParent(ofNode);
  while(thisNode && (ofNode==thisNode->left))
    {
      ofNode=thisNode;
      thisNode=getParent(thisNode);
    }
  return thisNode;
}



template<class T>
TigrBaseBinTreeNode<T> *TigrBinaryTree<T>::subTreeMax(TigrBaseBinTreeNode<T> *
						      subTree) const
{
  if(!subTree) 
    return NULL;

  while(subTree->right) 
    subTree=subTree->right;

  return subTree;
}



template<class T>
TigrBaseBinTreeNode<T> *TigrBinaryTree<T>::subTreeMin(TigrBaseBinTreeNode<T> *
						      subTree) const
{
  if(!subTree) 
    return NULL;

  while(subTree->left) 
    subTree=subTree->left;

  return subTree;
}



template<class T>
TigrBaseBinTreeNode<T> *TigrBinaryTree<T>::successor(TigrBaseBinTreeNode<T> *
						     ofNode)
{
  TigrBaseBinTreeNode<T> *thisNode;
   
  if(ofNode->right) 
    return subTreeMin(ofNode->right);
   
  thisNode=getParent(ofNode);
  while(thisNode && (ofNode==thisNode->right))
    {
      ofNode=thisNode;
      thisNode=getParent(thisNode);
    }
  return thisNode;
}



template<class T>
bool TigrBinaryTree<T>::insert(T t)
{
  TigrBaseBinTreeNode<T> *n=new TigrBaseBinTreeNode<T>;
  n->element=t;
  if(insert(n)) return true;
  delete n;
  return false;
}



template<class T>
bool TigrBinaryTree<T>::insert(TigrBaseBinTreeNode<T> *newNode)
{
  if(!root)
    {	
      root=newNode;
      newNode->left=newNode->right=NULL;
      cardinality=1;
      return true;
    }
   
  TigrBaseBinTreeNode<T> *thisNode=root;
  while(true) 
    {
      if(comparator.less(newNode->element,thisNode->element))
	{
	  if(!thisNode->left)
	    {
	      thisNode->left=newNode;
	      newNode->left=newNode->right=NULL;
	      ++cardinality;
	      return true;
	    }
	  thisNode=thisNode->left; 
	  continue;
	}
      else if(comparator.equal(newNode->element,thisNode->element))
	return false;
      else
	{
	  if(!thisNode->right)
	    {
	      thisNode->right=newNode;
	      newNode->left=newNode->right=NULL;
	      ++cardinality;
	      return true;
	    }
	  thisNode=thisNode->right; 
	  continue;
	}
    }
}



template<class T>
unsigned TigrBinaryTree<T>::getCardinality()
{
  return cardinality;
}



template<class T>
void TigrBinaryTree<T>::deleteNode(TigrBaseBinTreeNode<T> *thisNode)
{
  int numChildren;
  TigrBaseBinTreeNode<T> *parent, *child, *mySuccessor, *hisParent;
  numChildren=0;
  if(thisNode->left) numChildren++;
  if(thisNode->right) numChildren++;
  parent=getParent(thisNode);

  switch(numChildren) 
    {
    case 0:
      if(parent)
	if(parent->left==thisNode)
	  parent->left=NULL;
	else
	  parent->right=NULL;
      else  
	root=NULL;
      break;
      
    case 1:
      if(thisNode->left)
	{
	  child=thisNode->left;
	  thisNode->left=NULL;
	}
      else
	{
	  child=thisNode->right;
	  thisNode->right=NULL;
	}
      if(parent)
	if(parent->left==thisNode)
	  parent->left=child;
	else
	  parent->right=child;
      else  
	root=child;
      break;

    case 2:
      mySuccessor=thisNode->right;
      hisParent=NULL;
      while(mySuccessor->left)
	{
	  hisParent=mySuccessor;
	  mySuccessor=mySuccessor->left;
	}
      
      if(hisParent)  
	{
	  hisParent->left=mySuccessor->right;
	  mySuccessor->right=thisNode->right;
	}

      mySuccessor->left=thisNode->left;
      if(parent)
	if(parent->left==thisNode)
	  parent->left=mySuccessor;
	else
	  parent->right=mySuccessor;
      else  
	root=mySuccessor;
      thisNode->left=thisNode->right=NULL;
      break;
    }

  delete thisNode;
  --cardinality;
}



template<class T>
void TigrBinaryTree<T>::purge()
{
  if(!root) return;

  TigrBinTreeIter<T> iterator(this,TigrBinTreeIter<T>::POSTORDER);
  while(true)
    {
      if(!iterator) break;
      cout<<"deleting "<<*iterator<<endl;
      delete iterator.thisNode;
      ++iterator;
    }

  cardinality=0;
  root=NULL;
}



template<class T>
TigrBinaryTree<T>::iterator TigrBinaryTree<T>::begin()
{
  return TigrBinTreeIter(this);
}



template<class T>
TigrBinaryTree<T>::iterator TigrBinaryTree<T>::end()
{
  return TigrBinTreeIter(this,NULL,NULL,TigrBinTreeIter::FROM_NOWHERE,
			 TigrBinTreeIter::INORDER);
}



template<class T>
bool TigrBinaryTree<T>::isEmpty() const
{
  return cardinality==0;
}



template<class T>
TigrBaseBinTreeNode<T> *TigrBinaryTree<T>::find(T t)
{
  TigrBaseBinTreeNode<T> *n=root;
  while(n)
    if(comparator.less(t,n->element)) n=n->left;
    else if(comparator.equal(t,n->element)) return n;
    else n=n->right;
  return NULL;
}



template<class T>
void TigrBinaryTree<T>::remove(T t)
{
  TigrBaseBinTreeNode<T> *victim=find(t);
  if(victim) deleteNode(victim);
}



template<class T>
bool TigrBinaryTree<T>::getPredecessor(T of,T &is)
{
  TigrBaseBinTreeNode<T> *tp=find(of);
  if(!tp) throw "Element not found in TigrBinaryTree<T>::predecessor(T)";
  tp=predecessor(tp);
  if(!tp) return false;
  is=tp->element;
  return true;
}



template<class T>
bool TigrBinaryTree<T>::getSuccessor(T of,T &is)
{
  TigrBaseBinTreeNode<T> *tp=find(of);
  if(!tp) throw "Element not found in TigrBinaryTree<T>::successor(T)";
  tp=successor(tp);
  if(!tp) return false;
  is=tp->element;
  return true;
}


//==========================================================================
//                        TigrBinTreeIter methods
//==========================================================================



template<class T>
TigrBinaryTree<T>::TigrBinTreeIter::TigrBinTreeIter(TigrBinaryTree<T> *ofTree,
						    IterationType iterationType)
  : iterationType(iterationType),
    myTree(ofTree)
{
  reset();

  switch(iterationType)
    {
    case INORDER: thisNode=inOrder(); break;
    case PREORDER: thisNode=preOrder(); break;
    case POSTORDER: thisNode=postOrder(); break;
    }
}



template <class T>
TigrBinaryTree<T>::TigrBinTreeIter::TigrBinTreeIter(TigrBinaryTree<T> *tree,
				   TigrBaseBinTreeNode<T> *thisNode,
				   TigrBaseBinTreeNode<T> *nextNode,
				   TigrBinTreeIter::FromType fromWhere,
				   TigrBinTreeIter::IterationType it)
  : iterationType(it),
    myTree(tree),
    thisNode(thisNode),
    nextNode(nextNode),
    fromWhere(fromWhere)
{
}



template<class T>
TigrBaseBinTreeNode<T> *TigrBinaryTree<T>::TigrBinTreeIter::inOrder()
{
  TigrBaseBinTreeNode<T> *returnValue;
  
  while(true) 
    {
      if(nextNode==NULL) return NULL;
      switch(fromWhere) 
	{
	case FROM_LEFT:
	  returnValue=nextNode;
	  if(nextNode->right)
	    {
	      nextNode=nextNode->right;
	      fromWhere=FROM_PARENT;
	    }
	  else
	    {
	      nextNode=myTree->getParent(nextNode);
	      if(nextNode!=NULL)
		if(nextNode->left!=returnValue)
		  fromWhere=FROM_RIGHT;
	    }
	  return returnValue;
	  
	case FROM_RIGHT:
	  if(!myTree->getParent(nextNode)) 
	    {
	      nextNode=NULL;
	      return NULL;
	    }
	  if(myTree->getParent(nextNode)->left==nextNode)
	    fromWhere=FROM_LEFT;
	  nextNode=myTree->getParent(nextNode);
	  break;
	  
	case FROM_PARENT:
	  if(nextNode->left)
	    nextNode=nextNode->left;
	  else
	    fromWhere=FROM_LEFT;
	  break;
	  
	case FROM_NOWHERE:
	  return NULL;
	}
    }
}



template<class T>
TigrBaseBinTreeNode<T> *TigrBinaryTree<T>::TigrBinTreeIter::postOrder()
{
  TigrBaseBinTreeNode<T> *returnValue;

  while(true) 
    {
      if(nextNode==NULL) return NULL;
      switch(fromWhere) 
	{
	case FROM_RIGHT:
	  returnValue=nextNode;
	  nextNode=myTree->getParent(nextNode);
	  if(nextNode!=NULL)
	    if(nextNode->left==returnValue)
	      fromWhere=FROM_LEFT;
	  return returnValue;
	  
	case FROM_LEFT:
	  if(nextNode->right)
	    {
	      nextNode=nextNode->right;
	      fromWhere=FROM_PARENT;
	    }
	  else
	    fromWhere=FROM_RIGHT;
	  break;
	  
	case FROM_PARENT:
	  if(nextNode->left)
	    nextNode=nextNode->left;
	  else
	    fromWhere=FROM_LEFT;
	  break;
	  
	case FROM_NOWHERE:
	  return NULL;
	}
    }
}



template<class T>
TigrBaseBinTreeNode<T> *TigrBinaryTree<T>::TigrBinTreeIter::preOrder()
{
  TigrBaseBinTreeNode<T> *returnValue;

  while(true) 
    {
      if(nextNode==NULL) return NULL;
      switch(fromWhere) 
	{
	case FROM_PARENT:
	  returnValue=nextNode;
	  if(nextNode->left)
	    nextNode= nextNode->left;
	  else 
	    fromWhere=FROM_LEFT;
	  return returnValue;
	  
	case FROM_LEFT:
	  if(nextNode->right)
	    {
	      nextNode=nextNode->right;
	      fromWhere=FROM_PARENT;
	    }
	  else
	    fromWhere=FROM_RIGHT;
	  break;
	  
	case FROM_RIGHT:
	  returnValue=nextNode;
	  nextNode=myTree->getParent(nextNode);
	  if(nextNode!=NULL)
	    if(nextNode->left==returnValue)
	      fromWhere=FROM_LEFT;
	    else
	      fromWhere=FROM_RIGHT;
	  break;
	  
	case FROM_NOWHERE:
	  return NULL;
	}
    }
}



template<class T>
void TigrBinaryTree<T>::TigrBinTreeIter::reset()
{
  nextNode=myTree->getRoot();
  fromWhere=FROM_PARENT;
}



template<class T>
T &TigrBinaryTree<T>::TigrBinTreeIter::operator*()
{
  return thisNode->element;
}



template<class T>
TigrBinaryTree<T>::TigrBinTreeIter &TigrBinaryTree<T>::TigrBinTreeIter::
  operator++()
{
  switch(iterationType)
    {
    case PREORDER:   thisNode=preOrder();  break;
    case POSTORDER:  thisNode=postOrder(); break;
    case INORDER:    thisNode=inOrder();   break;
    }
  return *this;
}



template<class T>
TigrBinaryTree<T>::TigrBinTreeIter::operator bool()
{
  return thisNode!=NULL;
}



template<class T>
bool TigrBinaryTree<T>::TigrBinTreeIter::operator==(const 
						    TigrBinTreeIter<T> &other)
  const
{
  return thisNode==other.thisNode;
}



#endif
