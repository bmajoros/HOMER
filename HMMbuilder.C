/**************************************************************
HMMbuilder.C
bmajoros@tigr.org

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/

#include "HMMbuilder.H"

HiddenMarkovModel *HMMbuilder::randomHMM(int numStates,
					 float transitionDensity,
					 Alphabet &alphabet)
{
  HiddenMarkovModel &hmm=
    *new HiddenMarkovModel(alphabet,numStates,0);
  for(int i=0 ; i<numStates-1 ; ++i) hmm.addTransition(i,i+1);
  hmm.addTransition(numStates-1,0);
  for(int i=0 ; i<numStates ; ++i)
    for(int j=0 ; j<numStates ; ++j)
      if(Random0to1()<transitionDensity)
	hmm.addTransition(i,j);
  hmm.normalizeTransitions();
  return &hmm;
}

