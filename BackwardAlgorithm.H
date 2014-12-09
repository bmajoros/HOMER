/****************************************************************
 BackwardAlgorithm.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_BackwardAlgorithm_H
#define INCL_BackwardAlgorithm_H
#include "HiddenMarkovModel.H"
#include "tigr++/Sequence.H"
#include "tigr++/TigrArray1D.H"
#include "tigr++/TigrArray2D.H"

/**********************************************************************
 This is "Backward Algorithm."  It computes 
    b(k,i) = P[an HMM M will next emit x(i+1)..x(L) and then terminate
               | M is currently in state k]
 Compare this to the probability computed by the Forward Algorithm.
 The HMM begins in state 0 (which emits nothing) and continues its
 computation until it enters zero again (which again emits nothing).
 *********************************************************************/

class BackwardAlgorithm
{
public:
  BackwardAlgorithm(const TigrArray1D<double> &scalingFactors,
		    HiddenMarkovModel &,const TigrArray2D<int> &ngramArray,
		    int sequenceLength);
  double operator()(int state,int sequenceIndex); // Uses 1-based indices
  const TigrArray1D<double> &getScalingFactors();
  double getScaledP();
private:
  const int numStates, seqLen;
  TigrArray1D<double> scalingFactors;
  TigrArray2D<double> dpMatrix;
  double P, sum;
  HiddenMarkovModel &hmm;
  const TigrArray2D<int> &ngramArray;

  void compute();
};


#endif
