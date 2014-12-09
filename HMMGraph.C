/**************************************************************
HMMGraph.C
bmajoros@tigr.org

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/

#include "HMMGraph.H"
#include <iostream>
#include <math.h>
using namespace std;


HMMGraph::HMMGraph(HiddenMarkovModel &hmm)
  : hmm(hmm), predecessors(hmm.countStates()),
    emissions(hmm.getAlphabet().getNumElements())
{
  // ctor

  init();
}



void HMMGraph::init()
{
  int numStates=hmm.countStates();

  for(int i=0 ; i<numStates ; ++i)
    for(int j=0 ; j<numStates ; ++j)
      {
	double p=hmm.getTransitionProb(i,j);
	if(p>0)
	  predecessors[j].push_back(StateDoublePair(i,log(p),p));
      }

  Alphabet &alphabet=hmm.getAlphabet();
  int numSymbols=alphabet.getNumElements();
  for(int state=0 ; state<numStates ; ++state)
    for(int symbol=0 ; symbol<numSymbols ; ++symbol)
      {
	double p=hmm.getEmissionProb(state,symbol);
	if(p>0)
	  emissions[symbol].push_back(StateDoublePair(state,log(p),p));
      }
}



TigrVector<StateDoublePair> &HMMGraph::statesEmitting(Symbol s)
{
  return emissions[s];
}



TigrVector<StateDoublePair> &HMMGraph::statesPreceding(int state)
{
  return predecessors[state];
}



HiddenMarkovModel &HMMGraph::getHMM()
{
  return hmm;
}

