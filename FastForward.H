/**************************************************************
FastForward.H
bmajoros@tigr.org 1/1/2003

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_FastForward_H
#define INCL_FastForward_H
#include "HMMGraph.H"
#include "tigr++/TigrArray2D.H"
#include "tigr++/TigrArray1D.H"
#include "HiddenMarkovModel.H"
#include "tigr++/Sequence.H"

class HiddenMarkovModel;

class FastForward
{
public:
  FastForward(HiddenMarkovModel &,HMMGraph &,Sequence &);
  FastForward(HiddenMarkovModel &,HMMGraph &,Sequence &,int begin,int len);
  double getLogP();
  const TigrDblArray1D &getScalingFactors() const;
  bool PValueExceeds(FastForward &);
  double operator()(int state,int pos); // uses 1-based indices...
  double getScaledP();
private:
  double P;
  Sequence &sequence;
  int numStates, seqLen;
  TigrDblArray1D scalingFactors;
  TigrArray2D<double> dpMatrix;
  HiddenMarkovModel &hmm;
  HMMGraph &hmmGraph;
  const int offset;
  void computeDPMatrix();
  double computeScalingFactor(int);
};


#endif
