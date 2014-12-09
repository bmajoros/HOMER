/**************************************************************
HMMbuilder.H
bmajoros@tigr.org

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_HMMbuilder_H
#define INCL_HMMbuilder_H

#include "HiddenMarkovModel.H"
#include <iostream>
#include "tigr++/TigrRandom.H"
using namespace std;


class HMMbuilder
{
public:
  HiddenMarkovModel *randomHMM(int numStates,float transitionDensity,
			       Alphabet &);
};


#endif
