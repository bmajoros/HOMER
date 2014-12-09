/**************************************************************
HMMGraph.H
bmajoros@tigr.org

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/

#ifndef INCL_HMMGraph_H
#define INCL_HMMGraph_H

#include "HiddenMarkovModel.H"
#include "tigr++/TigrVector.H"

struct StateDoublePair
{
  short state;
  double logP, P;
  StateDoublePair(int s,double logP,double P) 
    : state(s), logP(logP), P(P) {}
};

class HMMGraph
{
public:
  HMMGraph(HiddenMarkovModel &);
  TigrVector<StateDoublePair> &statesEmitting(Symbol);
  TigrVector<StateDoublePair> &statesPreceding(int state);
  HiddenMarkovModel &getHMM();
private:
  HiddenMarkovModel &hmm;
  typedef TigrVector< TigrVector<StateDoublePair> > StateLists;
  StateLists emissions; // symbol -> states emitting it
  StateLists predecessors; // state -> predecessors
  void init();
};


#endif
