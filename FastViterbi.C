/**************************************************************
FastViterbi.C
bmajoros@tigr.org
***************************************************************/
#include <math.h>
#include <float.h>
#include <iostream>
#include "FastViterbi.H"
#include "tigr++/Constants.H"
using namespace std;

#define MOST_NEGATIVE -DBL_MAX



inline double safeAdd(double a,double b)
{
  if(a==MOST_NEGATIVE || b==MOST_NEGATIVE) return MOST_NEGATIVE;
  else return a+b;
}



FastViterbi::FastViterbi(HiddenMarkovModel &hmm)
  : numStates(hmm.countStates()), hmmGraph(hmm)
{
  // ctor
}



/**********************************************************************
  This method reconstructs the optimal path by performing traceback.
 */
TigrVector<int> *FastViterbi::buildPath(int L,int finalState,
					TigrArray2D<short> &linkBack)
{
  int currentState=finalState;
  TigrVector<int> &path=*new TigrVector<int>(L);
  path[L-1]=finalState;
  for(int pos=L-1 ; pos>0 ; --pos)
    {
      currentState=linkBack[currentState][pos+1];
      path[pos-1]=currentState;
    }
  return &path;
}



/**********************************************************************
  This method runs the dynamic programming algorithm by filling out the
  DP matrix and choosing the starting point for the optimal path, but
  it does not reconstruct the path (the caller can do that if descired).
 */
void FastViterbi::dp(const TigrArray2D<int> &ngramArray,
		     int &bestState,
		     double &bestScore,
		     TigrArray2D<double> &m,
		     TigrArray2D<short> &linkBack)
{
  // Fill out dynamic-programming array
  int L=ngramArray.getFirstDim();
  int Lminus1=L-1;
  fillOutArray(L,ngramArray,m,linkBack);

  // Choose best ending state
  bestState=1;
  bestScore=NEGATIVE_INFINITY;
  for(int i=0 ; i<numStates ; ++i)
    if(m[i][Lminus1]>bestScore)
      {
	bestState=i;
	bestScore=m[i][Lminus1];
      }
}



/**********************************************************************
  This method runs the dynamic programming algorithm and returns the
  score of the best path, but does not explicitly construct the best
  path (to save time).
 */
double FastViterbi::bestPathScore(const TigrArray2D<int> &ngramArray)
{
  int bestFinalState, L=ngramArray.getFirstDim();
  TigrArray2D<double> m(numStates,L+1);
  TigrArray2D<short> linkBack(numStates,L+1);
  double bestScore;
  dp(ngramArray,bestFinalState,bestScore,m,linkBack);
  return bestScore;
}



/**********************************************************************
  This method runs the dynamic programming algorithm and returns the 
  state sequence of the most probable path (not including the start 
  state).
 */
TigrVector<int> *FastViterbi::getPath(const TigrArray2D<int> &ngramArray,
				      double &pathScore)
{
  // Run the dynamic programming algorithm
  int bestFinalState, L=ngramArray.getFirstDim();
  int Lplus1=L+1;
  TigrArray2D<double> m(numStates,Lplus1);
  TigrArray2D<short> linkBack(numStates,Lplus1);
  dp(ngramArray,bestFinalState,pathScore,m,linkBack);
  
  // Reconstruct path via traceback
  return buildPath(L,bestFinalState,linkBack);
}



/**********************************************************************
  This method fills out the dynamic programming matrix, and that's all;
  it doesn't choose a best path (the caller can do that if necessary).
 */
void FastViterbi::fillOutArray(int L,
			       const TigrArray2D<int> &ngramArray,  
			       TigrArray2D<double> &m,
			       TigrArray2D<short> &linkBack)
{
  HiddenMarkovModel &hmm=hmmGraph.getHMM();
  int order=hmm.getOrder();
  int Lplus1=L+1;
  m.setAllTo(NEGATIVE_INFINITY);
  m[0][0]=0.0;
  for(int i=1 ; i<=L ; ++i)
    {
      Symbol s=ngramArray[i-1][0];//sequence[i-1];
      TigrVector<StateDoublePair> &statesEmittingS=
	hmmGraph.statesEmitting(s);
      int numEmitting=statesEmittingS.size();
      for(int l=0 ; l<numEmitting ; ++l)
	{
	  StateDoublePair &currentPair=statesEmittingS[l];
	  TigrVector<StateDoublePair> &precedingStates=
	    hmmGraph.statesPreceding(currentPair.state);
	  int numPreceding=precedingStates.size();
	  double bestP=log(0.0);
	  int bestPredecessor;
	  for(int j=0 ; j<numPreceding ; ++j)
	    {
	      StateDoublePair precedingPair=precedingStates[j];
	      double inductiveP=m[precedingPair.state][i-1];
	      double newP=safeAdd(inductiveP,precedingPair.logP);
	      if(newP>bestP)
		{
		  bestP=newP;
		  bestPredecessor=precedingPair.state;
		}
	    }

	  // Use the ngramArray-----------------------------------
	  double emissionProb=0;
	  for(int o=order ; o>=0 ; --o)
	    { // find the highest order with nonzero probability
	      int ngramId=ngramArray[i-1][o];
	      double p=hmm.getEmissionProb(currentPair.state,ngramId);
	      if(p>0) {emissionProb=p;break;}
	    }
	  m[currentPair.state][i]=safeAdd(bestP,log(emissionProb));

	  linkBack[currentPair.state][i]=bestPredecessor;
	}
    }

}
