/**************************************************************
HiddenMarkovModel.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_HiddenMarkovModel_H
#define INCL_HiddenMarkovModel_H
#include <iostream>
#include "tigr++/TigrString.H"
#include "tigr++/Sequence.H"
#include "tigr++/TigrArray2D.H"
#include "HigherOrderAlphabet.H"
using namespace std;


class HiddenMarkovModel
{
public:
  HiddenMarkovModel(const TigrString &fname,Alphabet &);
  HiddenMarkovModel(istream &,Alphabet &);
  HiddenMarkovModel(Alphabet &,int numStates,int order);
  HiddenMarkovModel(const HiddenMarkovModel &);
  double getEmissionProb(int inState,int symbol);
  double highestOrderEmissionProb(int inState,const TigrString &sequence,
				  int pos,HigherOrderAlphabet &);
  double highestOrderEmissionProb(int inState,
				  const TigrArray2D<int> &ngramArray,
				  int pos);
  void setTransitionProb(int fromState,int toState,double p);
  void setEmissionProb(int inState,int symbol,double p);
  void adjustEmissionProb(int inState,int symbol,double delta);
  bool doesTransitionExist(int from,int to);
  int countStates();
  Alphabet &getAlphabet();
  void printOn(ostream &);
  double getTransitionProb(int fromState,int toState);
  virtual bool save(const TigrString &filename);
  virtual void load(const TigrString &filename);
  virtual bool save(ostream &);
  virtual void load(istream &);
  void reverseComp();
  void normalizeTransitions();
  void normalizeEmissions();
  void addTransition(int from,int to);
  double getLogP(Sequence &seq,int begin=0,int len=-1);
  int getOrder() {return order;}
  int getNumSymbols() {return numSymbols;}
  void useVariableOrder();
  const HiddenMarkovModel &operator=(const HiddenMarkovModel &);
private:
  TigrArray2D<double> transitionProb;
  int numStates, order, numSymbols;
  TigrArray2D<double> emissionProb;
  Alphabet &alphabet;
  bool variableOrder;

  void swap(Symbol,Symbol);
  void initEmissionProbs();
  int computeNumSymbols();
};

ostream &operator<<(ostream &,HiddenMarkovModel &);

#endif
