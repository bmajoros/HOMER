/**************************************************************
HMMreader.H
bmajoros@tigr.org

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_HMMreader_H
#define INCL_HMMreader_H

#include "tigr++/TigrString.H"
#include "tigr++/TigrRegex.H"
#include "tigr++/Alphabet.H"
#include "HiddenMarkovModel.H"
#include "HigherOrderAlphabet.H"

/***********************************************************
  0->1 : 0.3
  0->2 : 0.7
  1->1 : 0.9 
  1->3 : 0.1
  ...etc...
  state 1 : A=0.3 T=0.2 C=0.4 g=0.1
  state 2 : ...etc...

  Note that 0 is the start state AND the final state, and
  it does not emit symbols.  The "state" statements are
  optional.  Any state not having one will be initialized
  with uniform emission probabilities.
***********************************************************/

class HMMreader
{
  Alphabet &alphabet;
  HigherOrderAlphabet &hoa;

  struct Transition
  {
    int from, to;
    float P;
    Transition(int from,int to,float P) : from(from), to(to), P(P) {}
  };
  struct Emission
  {
    int state;
    char symbol;
    float P;
    Emission(int s,char sym,float P) : state(s), symbol(sym), P(P) {}
  };
  TigrRegex transRegex, emitRegex, assignRegex;

public:
  HMMreader(Alphabet &,HigherOrderAlphabet &);
  HiddenMarkovModel *read(const TigrString &filename,int order=0);
};


#endif
