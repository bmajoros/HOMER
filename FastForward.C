/**************************************************************
FastForward.C
bmajoros@tigr.org 1/1/2003

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#include "tigr++/Constants.H"
#include <math.h>
#include "tigr++/TigrString.H"
#include "FastForward.H"
#include <iostream>
#include "HiddenMarkovModel.H"
using namespace std;


FastForward::FastForward(HiddenMarkovModel &hmm,HMMGraph &hmmGraph,
			 Sequence &sequence,
			 int begin,int len)
  : scalingFactors(len+1),sequence(sequence),
    dpMatrix(hmm.countStates(),len+1),
    numStates(hmm.countStates()),hmm(hmm),seqLen(len),
    hmmGraph(hmmGraph),offset(begin)
{
  computeDPMatrix();
}



FastForward::FastForward(HiddenMarkovModel &hmm,HMMGraph &hmmGraph,
			 Sequence &sequence)
  : scalingFactors(sequence.getLength()+1),sequence(sequence),
    dpMatrix(hmm.countStates(),sequence.getLength()+1),
    numStates(hmm.countStates()),hmm(hmm),seqLen(sequence.getLength()),
    hmmGraph(hmmGraph), offset(0)
{
  computeDPMatrix();
}



double FastForward::getLogP()
{
  double sum=log(P);
  for(int i=1 ; i<=seqLen ; ++i) sum+=log(scalingFactors[i]);
  return sum;
}



double FastForward::getScaledP()
{
  return P;
}



const TigrDblArray1D &FastForward::getScalingFactors() const
{
  return scalingFactors;
}



double FastForward::operator()(int state,int pos)
{
  return dpMatrix[state][pos];
}



double FastForward::computeScalingFactor(int i)
{
  double scalingFactor=0, sum;
  for(int l=1 ; l<numStates ; ++l)
    {
      sum=0;
      for(int k=0 ; k<numStates ; ++k)
	sum+=hmm.getTransitionProb(k,l)*dpMatrix[k][i-1];
      scalingFactor+=sum*hmm.getEmissionProb(l,sequence[offset+i-1]);
    }
  if(isNaN(scalingFactor) || isInfinity(scalingFactor))
    throw TigrString("scaling value in forward algorithm = ")+scalingFactor;
  if(scalingFactor<0 || scalingFactor>1)
    cerr << "scaling value in forward algorithm = "<<scalingFactor<<endl;
  return scalingFactor;
}



bool FastForward::PValueExceeds(FastForward &f)
{
  double difference=log(P)-log(f.P);
  for(int i=1 ; i<=seqLen ; ++i)
    {
      double d=log(scalingFactors[i])-log(f.scalingFactors[i]);
      difference+=d;
    }
  if(isInfinity(difference)||isNaN(difference)) 
    cerr << "FastForward: difference=" << difference << endl;
  bool exceeds=(difference>0);
  return exceeds;
}



void FastForward::computeDPMatrix()
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
	  TigrVector<StateDoublePair> &preceding=
	    hmmGraph.statesPreceding(l);
	  sum=0;
	  TigrVector<StateDoublePair>::iterator cur=preceding.begin(), 
	    end=preceding.end();
	  for(; cur!=end ; ++cur)
	    {
	      int fromState=(*cur).state;
	      double P=(*cur).P;
	      sum+=P*dpMatrix[fromState][i-1];
	    }
	  double v=dpMatrix[l][i]=sum*(hmm.getEmissionProb(l,
				  sequence[offset+i-1])/s_sub_i) ;
	  if(v<0 || v>1 || isNaN(v) || isInfinity(v))
	    throw TigrString("cell (")+l+","+i+")="+v;
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


