/****************************************************************
 BackwardAlgorithm.C
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#include "BackwardAlgorithm.H"
#include <math.h>
#include <iostream>
#include "tigr++/TigrString.H"
using namespace std;


BackwardAlgorithm::BackwardAlgorithm(const TigrArray1D<double> 
				     &scalingFactors,
				     HiddenMarkovModel &hmm,
				     const TigrArray2D<int> &ngrams,
				     int seqLen)
  : scalingFactors(seqLen+1),
    hmm(hmm),
    dpMatrix(hmm.countStates(),seqLen+1),
    numStates(hmm.countStates()),
    seqLen(seqLen),
    ngramArray(ngrams)
{
  compute();
}



const TigrArray1D<double> &BackwardAlgorithm::getScalingFactors()
{
  return scalingFactors;
}



double BackwardAlgorithm::getScaledP()
{
  return P;
}



double BackwardAlgorithm::operator()(int s,int i)
{
  return dpMatrix[s][i];
}



void BackwardAlgorithm::compute()
{
  for(int k=0 ; k<numStates ; ++k) 
    dpMatrix[k][seqLen]=hmm.getTransitionProb(k,0);
  for(int i=seqLen-1 ; i>=0 ; --i)
    for(int k=0 ; k<numStates ; ++k)
      {
	sum=0.0;	
	for(int l=1 ; l<numStates ; ++l)
	  sum+=dpMatrix[l][i+1]
	    * hmm.highestOrderEmissionProb(l,ngramArray,i) //### i?
	    * hmm.getTransitionProb(k,l);
	double scalingValue=0.0;
	for(int kk=0 ; kk<numStates ; ++kk)
	  for(int l=1 ; l<numStates ; ++l)
	    scalingValue+=dpMatrix[l][i+1]
	      * hmm.highestOrderEmissionProb(l,ngramArray,i) //### i?
	      * hmm.getTransitionProb(kk,l);
	if(isnan(scalingValue) || isinf(scalingValue))
	  throw TigrString("Backard scaling value ")+(i+1)+"="+scalingValue;
	/*if(scalingValue<0 || scalingValue>1)
	  cerr << "Backward scaling value " << (i+1) << "=" 
	  << scalingValue << endl;*/
	scalingFactors[i+1]=scalingValue;
	double entry;
	if(scalingValue==0.0) entry=dpMatrix[k][i]=0.0;
	else entry=dpMatrix[k][i]=sum/scalingValue;
	if(isnan(entry) || isinf(entry) || entry<0 || entry>1)
	  throw TigrString("b(")+k+","+i+")="+entry;
      }
  P=dpMatrix[0][0];
}



