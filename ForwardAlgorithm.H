/**************************************************************
ForwardAlgorithm.H
bmajoros@tigr.org 1/1/2003

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_ForwardAlgorithm_H
#define INCL_ForwardAlgorithm_H
#include "tigr++/TigrArray2D.H"
#include "tigr++/TigrArray1D.H"
#include "HiddenMarkovModel.H"

class ForwardAlgorithm
{
public:
  ForwardAlgorithm(HiddenMarkovModel &,const TigrArray2D<int> &ngramArray,
		   int sequenceLength);
  double getLogP();
  TigrDblArray1D &getScalingFactors();
  bool PValueExceeds(ForwardAlgorithm &);
  double operator()(int state,int pos); // uses 1-based indices...
  double getScaledP();
private:
  double P;
  int numStates, seqLen;
  TigrArray1D<double> scalingFactors;
  TigrArray2D<double> dpMatrix;
  HiddenMarkovModel &hmm;
  const TigrArray2D<int> &ngramArray;

  void computeDPMatrix();
  double computeScalingFactor(int);
};


#endif
