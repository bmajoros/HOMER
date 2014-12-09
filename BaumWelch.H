/****************************************************************
 BaumWelch.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/

#ifndef INCL_BaumWelch_H
#define INCL_BaumWelch_H

#include "ForwardAlgorithm.H"
#include "BackwardAlgorithm.H"
#include "HigherOrderAlphabet.H"
#include <math.h>
#include <iostream>
#include <fstream>
#include <float.h>
#include "tigr++/TigrVector.H"
#include "tigr++/TigrArray1D.H"
#include "tigr++/Sequence.H"
#include "tigr++/TigrProgress.H"
using namespace std;


class BaumWelch
{
public:
  BaumWelch(HiddenMarkovModel &,long maxIterations,
	    TigrVector<Sequence*> &trainingSet,
	    TigrVector<TigrString> &trainingSequences,
	    int minSampleSize,int order,
	    HigherOrderAlphabet &,
	    bool useNormalized);
  HigherOrderAlphabet &getHigherOrderAlphabet();
private:
  Alphabet &alphabet;
  HigherOrderAlphabet &higherOrderAlphabet;
  TigrProgress progress;
  TigrArray1D<double> likelihoods, scalingRatios;
  TigrArray2D<double> A; // expected transition counts
  TigrArray2D<double> E; // expected emission counts
  HiddenMarkovModel &hmm;
  long maxIterations;
  double logLikelihoodThreshold, sum;
  TigrVector<Sequence*> &trainingSet;
  TigrVector<TigrString> &trainingSequences;
  int numHmmStates, numTrain, alphabetSize, order, minSampleSize;
  bool useNormalized;
  
  void reviseExpectedEmitCounts(int sequenceLength,
				ForwardAlgorithm &,BackwardAlgorithm &,
				const TigrArray2D<int> &ngramArray);
  void reviseExpectedTransCounts(int sequenceLength,ForwardAlgorithm &,
				 BackwardAlgorithm &,
				 const TigrArray2D<int> &ngramArray);
  void mainAlgorithm();
  void fillLowerOrders();
  void renormalize();
};


#endif
