/**************************************************************
FastViterbi.H
bmajoros@tigr.org
***************************************************************/
#ifndef INCL_FastViterbi_H
#define INCL_FastViterbi_H

#include <iostream>
#include "tigr++/TigrSet.H"
#include "HiddenMarkovModel.H"
#include "tigr++/TigrBitSet.H"
#include "HMMGraph.H"
#include "tigr++/TigrVector.H"
#include "tigr++/TigrArray2D.H"
using namespace std;


class FastViterbi
{
public:
  FastViterbi(HiddenMarkovModel &);
  TigrVector<int> *getPath(const TigrArray2D<int> &ngramArray,
			   double &pathScore);
  double bestPathScore(const TigrArray2D<int> &ngramArray);
  HMMGraph &getHmmGraph() {return hmmGraph;}
private:
  int numStates;
  HMMGraph hmmGraph;
  void fillOutArray(int L,
		    const TigrArray2D<int> &ngramArray,  
		    TigrArray2D<double> &m,
		    TigrArray2D<short> &linkBack);
  TigrVector<int> *buildPath(int L,int finalState,
			     TigrArray2D<short> &linkBack);
  void dp(const TigrArray2D<int> &ngramArray,
	  int &bestState,
	  double &bestScore,
	  TigrArray2D<double> &m,
	  TigrArray2D<short> &linkBack);
};


#endif
