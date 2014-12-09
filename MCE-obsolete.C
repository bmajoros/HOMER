


struct Undo
{
  HiddenMarkovModel *hmm;
  int state;
  int symbol;
  double p;
  Undo() {}
  Undo(HiddenMarkovModel *h,int st,int sym,double p)
    : hmm(h), state(st), symbol(sym), p(p) {}
};


  TigrVector<Undo> undoables;


  void selectNgram(int order,TigrString &ngram);
  bool perturb(const TigrString &ngram);
  void renormalize(HiddenMarkovModel &,int state,const TigrString &ngram,
		   int changedSymbol,double changedP);
  void getNgrams(TigrVector<NGRAM_ARRAY*> &ngramArrays,
		 TigrMap<int,float> &counts);
  void normalizeNgramCounts(TigrMap<int,float> &TP_ngrams,
			    TigrMap<int,float> &TN_ngrams,
			    TigrMap<int,float> &FP_ngrams,
			    TigrMap<int,float> &FN_ngrams);
  void undo();



void Application::selectNgram(int len,TigrString &ngram)
{
  throw "OBSOLETE";

  ngram="";
  for(int i=0 ; i<len ; ++i)
    {
      int base=RandomNumber(4);
      if(base==3) ++base; // don't bother with the N symbol....
      char c=alphabet.lookup(base);
      ngram+=c;
    }
}



bool Application::perturb(const TigrString &ngram)
{
  throw "OBSOLETE";

  // Choose 1 state in each HMM to perturb
  int state1=RandomNumber(numStates1);
  int state2=RandomNumber(numStates2);

  // Get the current emission probabilities for this ngram
  int symbol=hoa->lookup(ngram);
  double p1=hmm1->getEmissionProb(state1,symbol);
  double p2=hmm2->getEmissionProb(state2,symbol);
  if(p1==0 || p2==0) return false;

  // Choose direction & size of perturbation
  const double DELTA=0.1;
  double delta=DELTA;
  if(p1<p2) 
  //if(Random0to1()<0.5) 
    delta=-delta;

  // Install new emission probabilities
  //hmm1->setEmissionProb(state1,symbol,p1+delta);
  //hmm2->setEmissionProb(state2,symbol,p2+delta);

  // Re-normalize to restore sum-to-one constraints
  renormalize(*hmm1,state1,ngram,symbol,p1+delta);
  renormalize(*hmm2,state2,ngram,symbol,p2+delta);

  return true;
}



void Application::renormalize(HiddenMarkovModel &hmm,
			      int state,
			      const TigrString &ngram,
			      int changedSymbol,
			      double changedP)
{
  throw "OBSOLETE";

  TigrString history=ngram.substr(1,ngram.length()-1);
  double sum=0;
  for(int base=0 ; base<alphabetSize ; ++base)
    {
      char c=alphabet.lookup(base);
      TigrString ngram=history+c;
      int symbol=hoa->lookup(ngram);
      double p=hmm.getEmissionProb(state,symbol);
      undoables.push_back(Undo(&hmm,state,symbol,p));
      if(symbol==changedSymbol) p=changedP;
      sum+=p;
    }
  for(int base=0 ; base<alphabetSize ; ++base)
    {
      char c=alphabet.lookup(base);
      TigrString ngram=history+c;
      int symbol=hoa->lookup(ngram);
      double p=
	(symbol==changedSymbol ?
	 changedP :
	 hmm.getEmissionProb(state,symbol));
      hmm.setEmissionProb(state,symbol,p/sum);
    }
}


void Application::normalizeNgramCounts(TigrMap<int,float> &TP_ngrams,
				       TigrMap<int,float> &TN_ngrams,
				       TigrMap<int,float> &FP_ngrams,
				       TigrMap<int,float> &FN_ngrams)
{
  throw "OBSOLETE";

  int numNgrams=hoa.size();
  for(int ngram=0 ; ngram<numNgrams ; ++ngram)
    {
      float &TP=TP_ngrams[ngram];
      float &TN=TN_ngrams[ngram];
      float &FP=FP_ngrams[ngram];
      float &FN=FN_ngrams[ngram];
      float total=TP+TN+FP+FN;
      TP/=total;
      TN/=total;
      FP/=total;
      FN/=total;
    }
}



void Application::getNgrams(TigrVector<NGRAM_ARRAY*> &ngramArrays,
			    TigrMap<int,float> &counts)
{
  throw "OBSOLETE";

  TigrVector<NGRAM_ARRAY*>::iterator cur=ngramArrays.begin(), 
    end=ngramArrays.end();
  for(; cur!=end ; ++cur)
    {
      NGRAM_ARRAY &ngramArray=**cur;
      int len=ngramArray.getFirstDim();
      for(int i=0 ; i<len ; ++i)
	{
	  int prevNgram=-1;
	  NGRAM_COLUMN &column=ngramArray[i];
	  for(int o=0 ; o<=order ; ++o)
	    {
	      int ngram=column[o];
	      if(ngram==prevNgram) break;
	      if(counts.isDefined(ngram)) ++counts[ngram];
	      else counts[ngram]=1;
	      prevNgram=ngram;
	    }
	}
    }
}





void Application::undo()
{
  throw "OBSOLETE";

  int n=undoables.size();
  for(int i=0 ; i<n ; ++i)
    {
      Undo &u=undoables[i];
      u.hmm->setEmissionProb(u.state,u.symbol,u.p);
    }
  undoables.clear();
}
