/**************************************************************
ForwardAlgorithm.C
bmajoros@tigr.org 1/1/2003

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#include "tigr++/Constants.H"
#include <math.h>
#include "tigr++/TigrString.H"
#include "ForwardAlgorithm.H"
#include <iostream>
using namespace std;


ForwardAlgorithm::ForwardAlgorithm(HiddenMarkovModel &hmm,
				   const TigrArray2D<int> &ngramArray,
				   int sequenceLength)
  : scalingFactors(sequenceLength+1),
    dpMatrix(hmm.countStates(),sequenceLength+1),
    numStates(hmm.countStates()),hmm(hmm),
    seqLen(sequenceLength),
    ngramArray(ngramArray)
{
  computeDPMatrix();
}



double ForwardAlgorithm::getLogP()
{
  double sum=log(P);
  for(int i=1 ; i<=seqLen ; ++i) sum+=log(scalingFactors[i]);
  return sum;
}



double ForwardAlgorithm::getScaledP()
{
  return P;
}



TigrDblArray1D &ForwardAlgorithm::getScalingFactors()
{
  return scalingFactors;
}



double ForwardAlgorithm::operator()(int state,int pos)
{
  return dpMatrix[state][pos];
}



double ForwardAlgorithm::computeScalingFactor(int i)
{
  double scalingFactor=0, sum;
  for(int l=1 ; l<numStates ; ++l)
    {
      sum=0;
      for(int k=0 ; k<numStates ; ++k)
	sum+=hmm.getTransitionProb(k,l)*dpMatrix[k][i-1];
      scalingFactor+=
	sum*hmm.highestOrderEmissionProb(l,ngramArray,i-1);
    }
  if(isNaN(scalingFactor) || isInfinity(scalingFactor))
    throw TigrString("scaling value in forward algorithm = ")+scalingFactor;
  /*if(scalingFactor<0 || scalingFactor>1)
    cerr << "scaling value in forward algorithm = "<<scalingFactor<<endl;*/
  return scalingFactor;
}



bool ForwardAlgorithm::PValueExceeds(ForwardAlgorithm &f)
{
  double difference=log(P)-log(f.P);
  for(int i=1 ; i<=seqLen ; ++i)
    {
      double d=log(scalingFactors[i])-log(f.scalingFactors[i]);
      difference+=d;
    }
  if(isInfinity(difference)||isNaN(difference)) 
    cerr << "ForwardAlgorithm: difference=" << difference << endl;
  bool exceeds=(difference>0);
  return exceeds;
}



void ForwardAlgorithm::computeDPMatrix()
{
  P=0;
  double s_sub_i, sum;
  dpMatrix[0][0]=1;
  for(int i=1 ; i<=seqLen ; ++i) dpMatrix[0][i]=0;
  for(int k=1 ; k<numStates ; ++k) dpMatrix[k][0]=0;
  for(int i=1 ; i<=seqLen ; ++i)
    {
      scalingFactors[i]=s_sub_i=computeScalingFactor(i);
      for(int l=1 ; l<numStates ; ++l)
	{
	  sum=0;	
	  for(int k=0 ; k<numStates ; ++k) 
	    sum+=hmm.getTransitionProb(k,l)*dpMatrix[k][i-1];
	  double v=dpMatrix[l][i]=
	    sum*(hmm.highestOrderEmissionProb(l,ngramArray,i-1)
		 /s_sub_i);
	  if(v<0 || v>1 || isNaN(v) || isInfinity(v))
	    throw TigrString("fw cell (")+l+","+i+")="+v;
	}
    }
  for(int k=1 ; k<numStates ; ++k) 
    P+=hmm.getTransitionProb(k,0)*dpMatrix[k][seqLen];
  if(P==0)
    for(int k=1 ; k<numStates ; ++k)
      {
	cout<<"dpMatrix[k][seqLen]="<<dpMatrix[k][seqLen]<<endl;
	cout<<"hmm.getTransitionProb(k,0)="
	    <<hmm.getTransitionProb(k,0)<<endl;
      }
}



