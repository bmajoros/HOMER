/**************************************************************
BaumWelch.C : the Baum-Welch Expectation Maximization algorithm
bmajoros@tigr.org

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/

#include "BaumWelch.H"
#include "tigr++/NthOrderStringIterator.H"

BaumWelch::BaumWelch(HiddenMarkovModel &hmm,long maxIterations,
		     TigrVector<Sequence*> &trainingSet,
		     TigrVector<TigrString> &trainingSequences,
		     int minSampleSize,int order,
		     HigherOrderAlphabet &hoa,
		     bool useNormalized)
  : hmm(hmm), 
    numHmmStates(hmm.countStates()),
    alphabet(hmm.getAlphabet()),
    numTrain(trainingSet.size()),
    likelihoods(trainingSet.size()),
    trainingSet(trainingSet),
    A(hmm.countStates(),hmm.countStates()),
    scalingRatios(1),
    maxIterations(maxIterations),
    E(0,0),
    minSampleSize(minSampleSize),
    trainingSequences(trainingSequences),
    order(order),
    higherOrderAlphabet(hoa),
    logLikelihoodThreshold(0.001),
    useNormalized(useNormalized)
{
  // First, run the canonical Baum-Welch algorithm, but with a larger
  // ("higher order") alphabet formed by grouping letters into ngrams
  alphabetSize=higherOrderAlphabet.size();
  E.resize(hmm.countStates(),alphabetSize);
  mainAlgorithm();

  // Pad the lower-order ngrams by summing over the higher-order ngrams
  // according to their lower-order suffixes
  fillLowerOrders();

  // Transform the unconditional probabilities of higher-order ngrams
  // into conditional probabilities of zeroth-order symbols, to produce
  // an HMM with higher-order emissions:
  renormalize();
}



/***************************************************************
  This is the main Baum-Welch Expectation Maximization algorithm
***************************************************************/
void BaumWelch::mainAlgorithm()
{
  double logLikelihood, deltaLikelihood=logLikelihoodThreshold+1, 
    oldLogLikelihood=0;

  /*************************************************************
    Various initialization steps
  *************************************************************/
  int numEmptyStrings=0;
  for(int j=0 ; j<numTrain ; ++j)
    {
      Sequence &sequence=*trainingSet[j];
      if(sequence.getLength()==0) ++numEmptyStrings;
    }
  bool allowEmptyStrings=hmm.doesTransitionExist(0,0);
  hmm.setTransitionProb(0,0,0);
  progress.start(maxIterations);

  /***********************************************************
    Main loop
  ***********************************************************/
  for(int i=0 ; i<maxIterations ; ++i)
    {
      A.setAllTo(0); // don't use pseudocounts--they cause likelihood drops!
      E.setAllTo(0); // (unless you re-do the math to account for them...)
      for(int symbol=0 ; symbol<alphabetSize ; ++symbol) E[0][symbol]=0;
      
      double likelihoodCorrection=0.0;
      for(int j=0 ; j<numTrain ; ++j)
	{
	  Sequence &sequence=*trainingSet[j];
	  TigrString &seqStr=trainingSequences[j];
	  int sequenceLength=sequence.getLength();
	  TigrArray2D<int> ngramArray(sequenceLength,order+1);
	  higherOrderAlphabet.initNgramArray(ngramArray,seqStr);
	  if(sequenceLength>0)
	    {
	      /***********************************************************
		Run the Forward algorithm to compute f(k,i)=P[a model M in
		the initial state will emit sequence x(1)..x(i) and reside
		in state k when emitting xi at time i]:
	       ***********************************************************/
	      ForwardAlgorithm f(hmm,ngramArray,sequenceLength);
	      const TigrArray1D<double> &scalingValues=
		f.getScalingFactors();
	      likelihoods[j]=f.getScaledP();
	      for(int i=1 ; i<=sequenceLength ; ++i)
		likelihoodCorrection+=log(scalingValues[i]);

	      /***********************************************************
		Run the Backward algorithm to compute b(k,i)=P[M will 
		next emit the string x(i+1)..x(L) and then terminate |
		M is currently in state k]:
	       ***********************************************************/
	      BackwardAlgorithm b(f.getScalingFactors(),hmm,
				  ngramArray,sequenceLength);

	      /***********************************************************
		Compute scaling ratios so the scaling values can be
		properly cancelled during the updating of the A and E
		counts:
	       ***********************************************************/
	      scalingRatios.resize(sequenceLength+1);
	      scalingRatios[sequenceLength]=
		1/f.getScalingFactors()[sequenceLength];
	      for(int i=sequenceLength-1 ; i>0 ; --i)
		{
		  double r=scalingRatios[i+1] *
		    b.getScalingFactors()[i+1]/f.getScalingFactors()[i];
		  scalingRatios[i]=r;
		  if(isinf(r) || isnan(r) || r==0) 
		    cerr << "WARNING: scaling ratio "<<i<<"="
			 <<scalingRatios[i]; 
		}

	      /***********************************************************
		Update the expected number of transitions and emissions
		used in generating this sequence:
	       ***********************************************************/
	      reviseExpectedTransCounts(sequenceLength,f,b,ngramArray);
	      reviseExpectedEmitCounts(sequenceLength,f,b,ngramArray);
	    }
	}

      /********************************************************************
	Recompute the a's and e's from the A's and E's (i.e., probabilities 
	from counts)
      ********************************************************************/
      for(int k=1 ; k<numHmmStates ; ++k)
	{
	  double sum=0;
	  for(int s=0 ; s<higherOrderAlphabet.size() ; ++s)
	    sum+=E[k][s];
	  if(sum==0)
	    for(int s=0 ; s<higherOrderAlphabet.size() ; ++s)
	      hmm.setEmissionProb(k,s,0);
	  else
	    for(int s=0 ; s<higherOrderAlphabet.size() ; ++s)
	      {
		double p=E[k][s]/sum;
		hmm.setEmissionProb(k,s,p);
	      }
	}
      for(int k=0 ; k<numHmmStates ; ++k)
	{
	  sum=0;
	  for(int ll=0 ; ll<numHmmStates ; ++ll) sum+=A[k][ll];
	  if(sum==0) for(int ll=0 ; ll<numHmmStates ; ++ll) 
	    hmm.setTransitionProb(k,ll,0);
	  else for(int l=0 ; l<numHmmStates ; ++l)
	    {
	      double p=A[k][l]/sum;
	      hmm.setTransitionProb(k,l,p);
	    }
	}
      if(useNormalized) {renormalize(); cerr<<"normalizing...";}

      /*********************************************************************
	Compute the log likelihood.  This is only necessary during 
        development and debugging, to ensure that recent changes to the 
        code have not changed the invariant that likelihood should increase
        monotonically (possibly except for tiny fluctuations due to 
        rounding in a digital computer):
      *********************************************************************/
      double logLikelihood=likelihoodCorrection;
      for(int j=0 ; j<numTrain ; ++j)
	if(trainingSet[j]->getLength()>0)
	  logLikelihood+=log(likelihoods[j]);
      deltaLikelihood=logLikelihood-oldLogLikelihood;
      oldLogLikelihood=logLikelihood;
      cerr << "log(likelihood)=" << logLikelihood;
      if(i>0)
	{
	  cerr << " " << progress.getProgress(i) << endl; 
	  if(deltaLikelihood<0)
	    cerr << "^--WARNING!  LIKELIHOOD DECREASED--^              >>> "
		 << deltaLikelihood << " <<<" << endl;
	}
      else cerr<<endl;
    }

  /*************************************************************************
   Allow start-state-self-transitions only if empty strings are to be 
   permitted
   ************************************************************************/
  if(allowEmptyStrings)
    {
      double p=double(numEmptyStrings)/double(numTrain);
      hmm.setTransitionProb(0,0,p);
      double q=1-p;
      for(int state=1 ; state<numHmmStates ; ++state)
	{
	  double rescaled=hmm.getTransitionProb(0,state)*q;
	  hmm.setTransitionProb(0,state,rescaled);
	}
    }
}



/*
  Update the expected emission counts for the HMM for a training sequence.
 */
void BaumWelch::reviseExpectedEmitCounts(int sequenceLength,
					 ForwardAlgorithm &f,
					 BackwardAlgorithm &b,
					 const TigrArray2D<int> &ngramArray)
{
  double seqProb=f.getScaledP();
  for(int k=1 ; k<numHmmStates ; ++k)
    for(int i=1 ; i<=sequenceLength ; ++i)
      {
	double delta=
	  f.getScalingFactors()[i] 
	  * scalingRatios[i]
	  * f(k,i) 
	  * b(k,i) 
	  / seqProb;
	E[k][ngramArray[i-1][order]]+=delta;
      }
}



/*
  Update the expected transition counts for the HMM for a training sequence.
 */
void BaumWelch::reviseExpectedTransCounts(int sequenceLength,
					  ForwardAlgorithm &f,
					  BackwardAlgorithm &b,
			     const TigrArray2D<int> &ngramArray)
{
  double seqProb=f.getScaledP();
  for(int k=0 ; k<numHmmStates ; ++k)
    for(int l=1 ; l<numHmmStates ; ++l)
      {
	sum=0.0;
	for(int i=0 ; i<sequenceLength ; ++i)
	  {
	    sum+=
	      f(k,i) 
	      * hmm.getTransitionProb(k,l)
	      * hmm.highestOrderEmissionProb(l,ngramArray,i) //### i?
	      * b(l,i+1) 
	      / seqProb
	      * scalingRatios[i+1];
	    if(isinf(sum) || isnan(sum)) throw TigrString("sum=")+sum+
		   " in BaumWelch::reviseExpectedTransCounts()";
	  }
	A[k][l]+=sum;
      }
  for(int k=1 ; k<numHmmStates ; ++k)
    A[k][0]+=f(k,sequenceLength) * hmm.getTransitionProb(k,0) / seqProb *
      scalingRatios[sequenceLength] * f.getScalingFactors()[sequenceLength];
  A[0][0]=0;
}



HigherOrderAlphabet &BaumWelch::getHigherOrderAlphabet()
{
  return higherOrderAlphabet;
}



void BaumWelch::fillLowerOrders()
{
  for(int o=order ; o>0 ; --o)
    for(int k=0 ; k<numHmmStates ; ++k)
      {
	NthOrderStringIterator iter(o+1,alphabet);
	while(!iter.done())
	  {
	    TigrString ngram=iter.getNextString();
	    int symbol=higherOrderAlphabet.lookup(ngram);
	    double p=hmm.getEmissionProb(k,symbol);
	    TigrString lowerOrderNgram=ngram.substr(1,o);
	    int lowerOrderSymbol=
	      higherOrderAlphabet.lookup(lowerOrderNgram);
	    hmm.adjustEmissionProb(k,lowerOrderSymbol,p);
	  }
      }
}



void BaumWelch::renormalize()
{
  int baseAlphabetSize=alphabet.size();
  for(int o=0 ; o<=order ; ++o)
    for(int k=0 ; k<numHmmStates ; ++k)
      {
	NthOrderStringIterator iter(o,alphabet);
	while(!iter.done())
	  {
	    TigrString history=iter.getNextString();
	    double sum=0.0;
	    for(Symbol a=0 ; a<baseAlphabetSize ; ++a)
	      {
		char c=alphabet.lookup(a);
		TigrString ngram=history+c;
		int symbol=higherOrderAlphabet.lookup(ngram);
		sum+=hmm.getEmissionProb(k,symbol);
	      }
	    if(sum>0.0)
	      for(Symbol a=0 ; a<baseAlphabetSize ; ++a)
		{
		  char c=alphabet.lookup(a);
		  TigrString ngram=history+c;
		  int symbol=higherOrderAlphabet.lookup(ngram);
		  hmm.setEmissionProb(k,symbol,
				      hmm.getEmissionProb(k,symbol)/sum);
		}
	  }
      }
}



