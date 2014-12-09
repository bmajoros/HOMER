/****************************************************************
 HOMER : Higher-Order Markovian Eukaryotic Recognizer
 A simple HMM for demonstrating HMM training & decoding algorithms

 (named in honor of the great Homer Simpson)

 bmajoros@tigr.org

 Copyright (c)2004, Bill Majoros.  This software is governed by  
 the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
#include "tigr++/TigrCommandLine.H"
#include "tigr++/TigrSet.H"
#include "tigr++/TigrRegex.H"
#include "tigr++/TigrVector.H"
#include "tigr++/TigrMap.H"
#include "FastViterbi.H"
#include "tigr++/TigrFastaReader.H"
#include "tigr++/TigrFile.H"
#include "tigr++/DnaAlphabet.H"
#include "tigr++/NthOrderStringIterator.H"
#include "HigherOrderAlphabet.H"
using namespace std;


class Exon
{
  int begin, end, transcriptId;
  TigrString axisId;
  char strand;
public:
  Exon(int begin,int end,int transcriptId,TigrString axisId,char strand)
    : begin(begin), end(end), transcriptId(transcriptId), axisId(axisId), 
      strand(strand) {}
  void printOn(ostream &);
  int getTranscriptId() const {return transcriptId;}
  char getStrand() const {return strand;}
};
ostream &operator<<(ostream &,Exon &);


class Application
{
  DnaAlphabet alphabet;
  HigherOrderAlphabet *higherOrderAlphabet;
  TigrSet<int> exonStates, intergenicStates;
  TigrRegex deflineRegex;
  vector<HiddenMarkovModel*> submodels;
  int numSymbols;
  vector<int> submodelOrigins; // state ids in composite HMM
  int numSubmodels;
  TigrVector<Exon*> exons;
  TigrMap<int,int> exonsPerGene;
  int order;

  void addToSet(TigrSet<int> &,int fromState,int numStates);
  void loadSubmodels(const TigrString &filename);
  HiddenMarkovModel *loadHMM(const TigrString &);
  void produceOutput();
  double getPathScore(vector<int> &path,Sequence &,HiddenMarkovModel &);
public:
  Application();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
  {
    try
      {
	Application app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



Application::Application()
  : deflineRegex(">\\s*(\\S+)")
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    TigrCommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=3) 
      throw string("homer [-p] <*.hmm> <submodels.txt> <*.fasta>\n");
    TigrString hmmFilename=cmd.arg(0);
    TigrString submodelsFilename=cmd.arg(1);
    TigrString fastaFilename=cmd.arg(2);

    time_t t=time(NULL);
    loadSubmodels(submodelsFilename);
    HiddenMarkovModel &hmm=*loadHMM(hmmFilename);
    FastViterbi viterbi(hmm);
    order=hmm.getOrder();
    higherOrderAlphabet=new HigherOrderAlphabet(alphabet,order);

    TigrFastaReader fastaReader(fastaFilename);
    TigrString defline, sequence;
    int transcriptId=0;
    while(fastaReader.nextSequence(defline,sequence))
      {
	if(!deflineRegex.search(defline))
	  throw TigrString("bad defline: ")+defline;
	TigrString axisId=deflineRegex[1];
	Sequence seq(sequence,alphabet);
	TigrArray2D<int> ngramArray(seq.getLength(),order+1);
	higherOrderAlphabet->initNgramArray(ngramArray,sequence);
	double pathScore;
	vector<int> &path=*viterbi.getPath(ngramArray,pathScore);
	int len=path.size();
	int exonBegin=-1;
	bool inExon=false, inIntergenic=false;
	for(int i=0 ; i<len ; ++i)
	  {
	    int state=path[i];
	    if(!inExon && exonStates.isMember(state))
	      {
		if(inIntergenic)
		  {
		    inIntergenic=false;
		    ++transcriptId;
		  }
		inExon=true;
		exonBegin=i;
	      }
	    else if(inExon && !exonStates.isMember(state))
	      {
		exons.push_back(new Exon(exonBegin,i,transcriptId,
					 axisId,'+'));
		inExon=false;
		if(intergenicStates.isMember(state)) inIntergenic=true;
		++exonsPerGene[transcriptId];
	      }
	    else if(!inIntergenic && intergenicStates.isMember(state))
	      inIntergenic=true;
	    else if(inIntergenic && !intergenicStates.isMember(state))
	      inIntergenic=false;
	  }
	cout << "#score="<<getPathScore(path,seq,hmm)<<endl;
	delete &path;
      }
    cout << "##gff-version 2\n"
	 << "##date " << ctime(&t)
	 << "##source-version homer 1.0 (Nov 2004)\n"
	 << "##contact bmajoros@tigr.org\n"
	 << "##type DNA\n";
    produceOutput();
    return 0;
  }


/*
void Application::loadSubmodels(const TigrString &filename)
{
  submodelOrigins.push_back(0);
  int nextFreeState=1;
  TigrFile file(filename);
  numSubmodels=0;
  while(!file.eof())
    {
      TigrString line=file.readLine();
      if(file.eof()) break;
      vector<TigrString> &fields=*line.getFields("= \t\n[]");
      int numFields=fields.size();
      if(numFields==1)
	throw TigrString("Syntax error in submodel list file: ")+line;
      if(numFields>0)
	{
	  int modelId=fields[0].asInt();
	  TigrString secondField=fields[1];
	  TigrString modelFilename=fields[numFields-1];
	  HiddenMarkovModel *submodel=
	    new HiddenMarkovModel(modelFilename,alphabet);
	  while(modelId>=submodels.size()) submodels.push_back(NULL);
	  submodels[modelId]=submodel;
	  int submodelStates=submodel->countStates()-1;
	  ++numSubmodels;
	  if(secondField=="exon")
	    addToSet(exonStates,nextFreeState,submodelStates);
	  else if(secondField=="intergenic")
	    addToSet(intergenicStates,nextFreeState,submodelStates);
	  submodelOrigins.push_back(nextFreeState);
	  nextFreeState+=submodelStates;
	}
      delete &fields;
    }
}
*/



void Application::loadSubmodels(const TigrString &filename)
{
  submodelOrigins.push_back(0);
  int nextFreeState=1;
  TigrFile file(filename);
  numSubmodels=0;
  while(!file.eof())
    {
      TigrString line=file.readLine();
      if(file.eof()) break;
      vector<TigrString> &fields=*line.getFields("= \t\n,");
      int numFields=fields.size();
      if(numFields==1)
	throw TigrString("Syntax error in submodel list file: ")+line;
      if(numFields>0)
	{
	  TigrString type=fields[0];
	  TigrSet<int> *theSet=NULL;
	  if(type=="intergenic") theSet=&intergenicStates;
	  else if(type=="exon") theSet=&exonStates;
	  else throw TigrString("Bad type in submodels file: ")+type;
	  for(int i=1 ; i<numFields ; ++i)
	    {
	      int state=fields[i].asInt();
	      if(state>0) theSet->insert(state);
	      //cout<<"inserting "<<state<<" as "<<type<<endl;
	    }
	}
      delete &fields;
    }
}



void Application::addToSet(TigrSet<int> &theSet,int fromState,int numStates)
{
  int end=fromState+numStates;
  for(int i=fromState ; i<end ; ++i)
    theSet.insert(i);
}



HiddenMarkovModel *Application::loadHMM(const TigrString &filename)
{
  HiddenMarkovModel *H=new HiddenMarkovModel(filename,alphabet);
  return H;
}



void Exon::printOn(ostream &os)
{
  os << axisId << "\thomer\t" << "exon" << "\t"
     << begin+1 << "\t" << end << "\t.\t" << strand << "\t"
     << "." <<"\ttransgrp="<< transcriptId << ";" << endl;
}



ostream &operator<<(ostream &os,Exon &exon)
{
  exon.printOn(os);
  return os;
}




void Application::produceOutput()
{
  int n=exons.size();
  for(int i=0 ; i<n ; ++i)
    {
      Exon *exon=exons[i];
      cout << *exon;
    }
}




double Application::getPathScore(vector<int> &path,Sequence &seq,
				 HiddenMarkovModel &hmm)
{
  double logP=0;
  int pathLength=path.size();
  for(int i=0 ; i<pathLength ; ++i)
    {
      int state=path[i];
      Symbol symbol=seq[i];
      logP+=log(hmm.getEmissionProb(state,symbol));
      if(i+1<pathLength) 
	logP+=log(hmm.getTransitionProb(state,path[i+1]));
    }
  return logP;
}





