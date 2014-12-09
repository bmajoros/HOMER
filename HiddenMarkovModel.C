/**************************************************************
HiddenMarkovModel.C
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#include <fstream>
#include "HiddenMarkovModel.H"
#include "tigr++/TigrRandom.H"
#include <iostream>
#include "FastForward.H"
#include "HMMGraph.H"
using namespace std;


HiddenMarkovModel::HiddenMarkovModel(Alphabet &alphabet,int numStates,
				     int order)
  : numStates(numStates), 
    emissionProb(0,0),
    transitionProb(numStates,numStates), 
    alphabet(alphabet),
    order(order),
    variableOrder(false)
{
  numSymbols=computeNumSymbols();
  emissionProb.resize(numStates,numSymbols);

  emissionProb.setAllTo(0.0);
  transitionProb.setAllTo(0.0);
}



HiddenMarkovModel::HiddenMarkovModel(const HiddenMarkovModel &other)
  : emissionProb(other.emissionProb), 
    numStates(other.numStates), 
    transitionProb(other.transitionProb), 
    alphabet(other.alphabet),
    order(other.order),
    variableOrder(other.variableOrder),
    numSymbols(other.numSymbols)
{
}



const HiddenMarkovModel &HiddenMarkovModel::operator=(
                                    const HiddenMarkovModel &other)
{
  emissionProb=other.emissionProb;
  numStates=other.numStates;
  transitionProb=other.transitionProb; 
  order=other.order;
  variableOrder=other.variableOrder;
  numSymbols=other.numSymbols;
  return other;
}



HiddenMarkovModel::HiddenMarkovModel(const TigrString &filename,
				     Alphabet &alphabet)
  : emissionProb(0,0), 
    transitionProb(0,0), 
    alphabet(alphabet),
    variableOrder(false)
{
  load(filename);
}



HiddenMarkovModel::HiddenMarkovModel(istream &is,Alphabet &alphabet)
  : transitionProb(0,0), 
    emissionProb(0,0), 
    alphabet(alphabet),
    variableOrder(false)
{
  load(is);
}



void HiddenMarkovModel::printOn(ostream &os)
{
  os << "Transitions:\n";
  for(int i=0 ; i<numStates ; ++i)
    {
      for(int j=0 ; j<numStates ; ++j)
	os << '\t' << transitionProb[i][j];
      os << '\n';
    }
  int numSymbols=alphabet.getNumElements();
  os << "\nEmissions:\n";
  for(int i=0 ; i<numStates ; ++i)
    {
      for(int j=0 ; j<numSymbols ; ++j)
	os << '\t' << emissionProb[i][j];
      os << '\n';
    }  
  os << endl;
}



ostream &operator<<(ostream &os,HiddenMarkovModel &model)
{
  model.printOn(os);
  return os;
}



double HiddenMarkovModel::getTransitionProb(int from,int to)
{
  return transitionProb[from][to];
}



void HiddenMarkovModel::load(const TigrString &fname)
{
  ifstream is(fname.c_str());
  if(!is.good()) throw TigrString("Can't open ")+fname;
  load(is);
}



int HiddenMarkovModel::computeNumSymbols()
{
  int a=alphabet.getNumElements(), n=0;
  for(int i=0 ; i<=order ; ++i)
    n+=(int) pow((double)a,(double)i+1);//###
  return n;
}



void HiddenMarkovModel::load(istream &is)
{
  numSymbols=alphabet.getNumElements();

  // Load number of states and allocate arrays of that size
  is >> numStates >> order;
  numSymbols=computeNumSymbols();
  emissionProb.resize(numStates,numSymbols);
  transitionProb.resize(numStates,numStates);

  // Read transition  && emission probabilities
  is >> transitionProb >> emissionProb;
}



void HiddenMarkovModel::normalizeTransitions()
{
  double sum;
  for(int i=0 ; i<numStates ; ++i)
    {
      sum=0;
      for(int j=0 ; j<numStates ; ++j)
	sum+=transitionProb[i][j];
      if(sum>0)
	for(int j=0 ; j<numStates ; ++j)
	  transitionProb[i][j]/=sum;
    }
}



bool HiddenMarkovModel::save(const TigrString &fname)
{
  ofstream os(fname.c_str());
  if(!os.good()) throw TigrString("Can't create ")+fname;
  bool success=save(os);
  return success;
}



bool HiddenMarkovModel::save(ostream &os)
{
  os << numStates << '\t' << order << '\n' 
     << transitionProb << '\n' 
     << emissionProb << endl;
  return true;
}



void HiddenMarkovModel::reverseComp()
{
  double t;
  Symbol C=alphabet.lookup('C');
  Symbol T=alphabet.lookup('T');
  Symbol A=alphabet.lookup('A');
  Symbol G=alphabet.lookup('G');
  swap(G,C);
  swap(A,T);
  normalizeEmissions();
  for(int i=0 ; i<numStates ; ++i)
    for(int j=i+1 ; j<numStates ; ++j)
      {
	t=transitionProb[j][i];
	transitionProb[j][i]=transitionProb[i][j];
	transitionProb[i][j]=t;
      }

  // NOTE:  do NOT call normalizeTransitions() here unless
  //        you really want to ... otherwise, the reverseComp
  //        operation will not be "reversible" (i.e., a
  //        symmetric binary relation).
}



bool HiddenMarkovModel::doesTransitionExist(int from,int to)
{
  bool exists=transitionProb[from][to]>0;
  return exists;
}



void HiddenMarkovModel::normalizeEmissions()
{
  throw "HiddenMarkovModel::normalizeEmissions()";
  double sum;
  int numSymbols=alphabet.getNumElements();
  for(int from=0 ; from<numStates ; ++from)
    {
      sum=0;
      for(int s=0 ; s<numSymbols ; ++s)
	sum+=emissionProb[from][s];
      if(sum>0)
	for(int s=0 ; s<numSymbols ; ++s)
	  emissionProb[from][s]/=sum;
    }
}



void HiddenMarkovModel::addTransition(int from,int to)
{
  setTransitionProb(from,to,Random0to1());
}



int HiddenMarkovModel::countStates()
{
  return numStates;
}



void HiddenMarkovModel::swap(Symbol s1,Symbol s2)
{
  double t;
  for(int state=0 ; state<numStates ; ++state)
    {
      t=emissionProb[state][s1];
      emissionProb[state][s1]=emissionProb[state][s2];
      emissionProb[state][s2]=t;
    }
}



void HiddenMarkovModel::setTransitionProb(int from,int to,double prob)
{
  transitionProb[from][to]=prob;
}



double HiddenMarkovModel::getEmissionProb(int state,int s)
{
  return emissionProb[state][s];
}



Alphabet &HiddenMarkovModel::getAlphabet()
{
  return alphabet;
}



double HiddenMarkovModel::getLogP(Sequence &seq,int begin,int len)
{
  if(len<0) len=seq.getLength();
  HMMGraph hmmGraph(*this);
  FastForward f(*this,hmmGraph,seq,begin,len);
  return f.getLogP();
}



void HiddenMarkovModel::initEmissionProbs()
{
  double sum;
  int numSymbols=alphabet.getNumElements();
  for(int i=0 ; i<numStates ; ++i)
    {
      sum=0.0;
      for(int j=0 ; j<numSymbols ; ++j)
	{
	  double r=RandomFloat(0.20,0.80);
	  sum+=r;
	  emissionProb[i][j]=r;
	}
      if(sum>0)
	for(int j=0 ; j<numSymbols ; ++j)
	  { emissionProb[i][j]/=sum; }
    }
}



void HiddenMarkovModel::setEmissionProb(int state,int s,double prob)
{
  emissionProb[state][s]=prob;
}



void HiddenMarkovModel::adjustEmissionProb(int inState,int symbol,
					   double delta)
{
  emissionProb[inState][symbol]+=delta;
}



double HiddenMarkovModel::highestOrderEmissionProb(int inState,
			    const TigrArray2D<int> &ngramArray,
						       int pos)
{
  const TigrArray2D<int>::RowIn2DArray<int> &row=ngramArray[pos];
  if(variableOrder)
    for(int i=order ; i>=0 ; --i)
      {
	int symbol=row[i];
	double p=getEmissionProb(inState,symbol);
	if(p>0.0) return p;
	// otherwise, P==0, so consider order i-1 in next iteration
      }
  else return getEmissionProb(inState,row[order]);
  return 0.0;
}



double HiddenMarkovModel::highestOrderEmissionProb(int inState,
				    const TigrString &sequence,
			   int pos,HigherOrderAlphabet &hAlpha)
{
  double p;
  if(variableOrder)
    {
      int begin=pos-order;
      if(begin<0) begin=0;
      for(; begin<=pos ; ++begin)
	{
	  int len=pos-begin+1;
	  TigrString ngram=sequence.substr(begin,len);
	  int symbol=hAlpha.lookup(ngram);
	  p=getEmissionProb(inState,symbol);
	  if(p>0.0) return p;
	  //else P==0, so consider a lower-order emission next iteration
	}
    }
  else
    {
      int begin=pos-order;
      if(begin<0) begin=0;
      int len=pos-begin+1;
      TigrString ngram=sequence.substr(begin,len);
      int symbol=hAlpha.lookup(ngram);
      p=getEmissionProb(inState,symbol);
    }
  return p;
}



void HiddenMarkovModel::useVariableOrder()
{
  variableOrder=true;
  cout<<"WARNING!  Using variable order HMM"<<endl;
}


