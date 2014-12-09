/**************************************************************
HMMreader.C
bmajoros@tigr.org

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#include "HMMreader.H"
#include "tigr++/TigrFile.H"
#include "tigr++/TigrRandom.H"
#include <iostream>
#include "HiddenMarkovModel.H"
using namespace std;


HMMreader::HMMreader(Alphabet &alphabet,HigherOrderAlphabet &hoa)
  : alphabet(alphabet),
    transRegex("(\\d+)\\s*->\\s*(\\d+)\\s*:\\s*(\\S+)"),
    emitRegex("state\\s*(\\d+)\\s*:\\s*(\\S.*\\S)"),
    assignRegex("(.)=(\\S+)"),
    hoa(hoa)
{
  // ctor
}



HiddenMarkovModel *HMMreader::read(const TigrString &filename,int order)
{
  TigrVector<Transition*> transitions;
  TigrVector<Emission*> emissions;

  int maxState=0;
  TigrFile file(filename);
  while(!file.eof())
    {
      TigrString line=file.readLine();
      line.trimWhitespace();
      if(transRegex.search(line))
	{
	  int from=transRegex[1].asInt(), to=transRegex[2].asInt();
	  float P=transRegex[3].asFloat();
	  if(from>maxState) maxState=from;
	  if(to>maxState) maxState=to;
	  transitions.push_back(new Transition(from,to,P));
	}
      else if(emitRegex.search(line))
	{
	  int state=emitRegex[1].asInt();
	  TigrString distrString=emitRegex[2];
	  TigrVector<TigrString> *fields=distrString.getFields();
	  TigrVector<TigrString>::iterator cur=fields->begin();
	  for( ; cur!=fields->end() ; ++cur)
	    {
	      if(assignRegex.search(*cur))
		{
		  char symbol=assignRegex[1][0];
		  float P=assignRegex[2].asFloat();
		  emissions.push_back(new Emission(state,symbol,P));
		}
	      else 
		{
		  (*cur).trimWhitespace();
		  if((*cur).length()>0)
		    throw TigrString("Syntax error in HMM structure file: ")
		      + *cur;
		}
	    }
	  delete fields;
	}
      else if(line.length()>0)
	throw TigrString("Syntax error in HMM structure file: ")+line;
    }
  
  HiddenMarkovModel *hmm=new HiddenMarkovModel(alphabet,maxState+1,order);
  TigrVector<Transition*>::iterator cur=transitions.begin();
  for(; cur!=transitions.end() ; ++cur)
    {
      Transition *trans=*cur;
      hmm->setTransitionProb(trans->from,trans->to,trans->P);
      delete trans;
    }

  if(emissions.size()>0)
    {
      if(order>0) throw "Error: HMM structure file contains emission probabilities, but HMM is not zeroth order\n";
      TigrVector<Emission*>::iterator ecur=emissions.begin();
      for(; ecur!=emissions.end() ; ++ecur)
	{
	  Emission *emission=*ecur;
	  hmm->setEmissionProb(emission->state,
			       alphabet.lookup(emission->symbol),
			       emission->P);
	  delete emission;
	}
    }
  else
    {
      /* 
      // UNIFORM:
      int numNgrams=hoa.size();
      double p=1.0/numNgrams;
      for(int s=0 ; s<numNgrams ; ++s)
	for(int i=0 ; i<=maxState ; ++i)
	  hmm->setEmissionProb(i,s,p);
      */

      // RANDOM:
      int numNgrams=hoa.size();
      double sum=0;
      for(int s=0 ; s<numNgrams ; ++s)
	for(int i=0 ; i<=maxState ; ++i)
	  {
	    double p=Random0to1();
	    hmm->setEmissionProb(i,s,p);
	    sum+=p;
	  }
      for(int i=0 ; i<=maxState ; ++i)
	for(int s=0 ; s<numNgrams ; ++s)
	  hmm->setEmissionProb(i,s,hmm->getEmissionProb(i,s)/sum);

    }
 
  return hmm;
}



