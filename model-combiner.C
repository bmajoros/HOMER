/**************************************************************
model-combiner.C
bmajoros@tigr.org

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/

#include <iostream>
#include <fstream>
#include "tigr++/TigrRandom.H"
#include "tigr++/TigrFile.H"
#include "tigr++/TigrCommandLine.H"
#include "tigr++/Alphabet.H"
#include "tigr++/Sequence.H"
#include "HiddenMarkovModel.H"
#include "HMMbuilder.H"
#include "HMMreader.H"
#include "tigr++/TigrFastaReader.H"
#include "tigr++/TigrMap.H"
using namespace std;


class Application
{
  Alphabet alphabet;
  HiddenMarkovModel *metaModel;
  int numSubmodels, numSymbols;
  TigrVector<Sequence*> *trainingSet;
  TigrMap<int,TigrMap<int,int> > trans;   // state x state x count
  TigrMap<int,TigrMap<Symbol,int> > emit; // state x symbol x count
  TigrVector<int> submodelOrigins; // state ids in composite HMM
  TigrVector<TigrString> submodelNames;
  
  TigrVector<HiddenMarkovModel*> *loadSubmodels(const TigrString &filename);
  void initAlphabet();
  HiddenMarkovModel *merge(TigrVector<HiddenMarkovModel*> &submodels);
  void copyIn(HiddenMarkovModel &submodel,HiddenMarkovModel &compositeHMM,
	      int startingAt);
  void linkModels(TigrVector<HiddenMarkovModel*> &,int fromModelId,
		  int toModelId,HiddenMarkovModel &compositeHMM,
		  double metaP);
  void closeTransitivity(int fromState,int toModelId,double fromP,
			 HiddenMarkovModel &compositeHMM,
			 TigrVector<HiddenMarkovModel*> &);
public:
  Application() {}
  int go(int argc,char *argv[]);
};


int main(int argc,char *argv[])
{
  try
    {
      Application app;
      return app.go(argc,argv);
    }
  catch(const char *msg)
    {
      cerr << msg << endl;
      return -1;
    }
  catch(const TigrString &s)
    {
      cerr << s.c_str() << endl;
      return -1;
    }
  catch(string s)
    {
      cerr << s.c_str() << endl;
      return -1;
    }
  catch(...)
    {
      cerr << "unknown exception caught in main()" << endl;
      return -1;
    }
  return 0;
}



void usage()
{
  throw "model-combiner <meta-model.hmms> <submodels.txt> <outfile.hmm>";
}


int Application::go(int argc,char *argv[])
{
  // Process command line
  TigrCommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3) usage();
  TigrString structureFile=cmd.arg(0);
  TigrString submodelListFile=cmd.arg(1);
  TigrString outfile=cmd.arg(2);
  randomize();

  // Load the MetaModel
  initAlphabet();
  numSymbols=alphabet.getNumElements();
  HMMreader reader(alphabet);
  metaModel=reader.read(structureFile);
  //cerr << "Metamodel:\n" << *metaModel << endl;

  // Load the submodels
  TigrVector<HiddenMarkovModel*> &submodels=
    *loadSubmodels(submodelListFile);
  numSubmodels=submodels.size()-1;
  submodels[0]=NULL;

  // Merge the submodels into a composite model
  HiddenMarkovModel *compositeHMM=merge(submodels);

  // Save the results
  compositeHMM->save(outfile);

  return 0;
}



TigrVector<HiddenMarkovModel*> *Application::loadSubmodels(
                                              const TigrString &filename)
{
  submodelNames.push_back("");
  TigrVector<HiddenMarkovModel*> *models=new TigrVector<HiddenMarkovModel*>;
  TigrFile file(filename);
  while(!file.eof())
    {
      TigrString line=file.readLine();
      if(file.eof()) break;
      TigrVector<TigrString> &fields=*line.getFields("= \t\n+,");
      int numFields=fields.size();
      if(numFields==1)
	throw TigrString("Syntax error in submodel list file: ")+line;
      if(numFields>0)
	{
	  TigrString firstField=fields[0];
	  int modelId=firstField.asInt();
	  TigrString modelFilename=fields[1];
	  delete &fields;
	  submodelNames.push_back(modelFilename);
	  
	  HiddenMarkovModel *submodel=
	    new HiddenMarkovModel(modelFilename,alphabet);
	  while(modelId>=models->size()) models->push_back(NULL);
	  (*models)[modelId]=submodel;
	}
    }
  return models;
}



void Application::initAlphabet()
{
  alphabet.add('A');
  alphabet.add('C');
  alphabet.add('G');
  alphabet.add('N');
  alphabet.add('T');
}



void Application::copyIn(HiddenMarkovModel &submodel,
			 HiddenMarkovModel &compositeHMM,int delta)
{
  // Copy transitions
  int numStates=submodel.countStates();
  for(int i=1 ; i<numStates ; ++i)
    {
      for(int j=1 ; j<numStates ; ++j)
	compositeHMM.setTransitionProb(i+delta,j+delta,
				  submodel.getTransitionProb(i,j));
      compositeHMM.setTransitionProb(i+delta,0,
				submodel.getTransitionProb(i,0));
    }
  
  // Copy emissions
  for(int i=1 ; i<numStates ; ++i)
    for(int j=0 ; j<numSymbols ; ++j)
      compositeHMM.setEmissionProb(i+delta,Symbol(j),
			       submodel.getEmissionProb(i,Symbol(j)));
}



HiddenMarkovModel *Application::merge(TigrVector<HiddenMarkovModel*> 
				      &submodels)
{
  // Determine how many states we need in the composite model
  int totalStates=1; // has a silent zero state
  //cout << numSubmodels << " submodels:" << endl;
  for(int i=1 ; i<=numSubmodels ; ++i)
    {
      HiddenMarkovModel &submodel=*submodels[i];
      totalStates+=submodel.countStates()-1;//don't count zero states
    }

  // Allocate a composite model
  HiddenMarkovModel &compositeHMM=
    *new HiddenMarkovModel(alphabet,totalStates);

  // Copy the submodels into distinct, unconnected regions of the
  // composite model
  int nextFreeState=1;
  submodelOrigins.push_back(0);
  for(int modelId=1 ; modelId<=numSubmodels ; ++modelId)
    {
      HiddenMarkovModel &submodel=*submodels[modelId];
      copyIn(submodel,compositeHMM,nextFreeState-1);
      submodelOrigins.push_back(nextFreeState);
      /*
      cout << "submodel #" << modelId << " has " 
	   << submodel.countStates()-1 
	   << " states, numbered " << nextFreeState << "-"
	   << nextFreeState+submodel.countStates()-2 
	   << "\t" << submodelNames[modelId] << endl;
      */
      nextFreeState+=submodel.countStates()-1;
    }
  /*
  cout << "total states=" << totalStates << ", alphabet=" 
       << alphabet << endl;
  */

  // Establish connections between the submodels, as dictated by
  // the MetaModel
  for(int from=1 ; from<=numSubmodels ; ++from)
    for(int to=1 ; to<=numSubmodels ; ++to)
      if(metaModel->doesTransitionExist(from,to))
	linkModels(submodels,from,to,compositeHMM,
		   metaModel->getTransitionProb(from,to));

  // Establish connections from state 0 in the metaModel
  for(int toModelId=1 ; toModelId<=numSubmodels ; ++toModelId)
    if(metaModel->doesTransitionExist(0,toModelId))
      {
	double metaP=metaModel->getTransitionProb(0,toModelId);
	HiddenMarkovModel &toModel=*submodels[toModelId];
	int toModelDelta=submodelOrigins[toModelId]-1;
	int numToModelStates=toModel.countStates();
	for(int toState=1 ; toState<numToModelStates ; ++toState)
	  {
	    if(toModel.doesTransitionExist(0,toState))
	      {
		double newP=metaP*toModel.getTransitionProb(0,toState);
		compositeHMM.setTransitionProb(0,toState+toModelDelta,newP);
	      }
	  }

	// Handle 0->0 transitions
	if(toModel.doesTransitionExist(0,0))
	  {
	    //cout << "MODEL "<<toModelId<<" HAS 0->0 TRANSITION"<<endl;
	    double newP=metaP*toModel.getTransitionProb(0,0);
	    closeTransitivity(0,toModelId,newP,compositeHMM,submodels);
	  }
      }

  // Establish connections to state 0 in the metaModel
  for(int fromModelId=1 ; fromModelId<=numSubmodels ; ++fromModelId)
    if(metaModel->doesTransitionExist(fromModelId,0))
      {
	double metaP=metaModel->getTransitionProb(fromModelId,0);
	HiddenMarkovModel &fromModel=*submodels[fromModelId];
	int fromModelDelta=submodelOrigins[fromModelId]-1;
	int numFromModelStates=fromModel.countStates();
	for(int fromState=1 ; fromState<numFromModelStates ; ++fromState)
	  {
	    if(fromModel.doesTransitionExist(fromState,0))
	      {
		double newP=metaP*fromModel.getTransitionProb(fromState,0);
		compositeHMM.setTransitionProb(fromState+fromModelDelta,0,
					       newP);
	      }
	  }
      }
    else
      {
	HiddenMarkovModel &fromModel=*submodels[fromModelId];
	int fromModelDelta=submodelOrigins[fromModelId]-1;
	int numFromModelStates=fromModel.countStates();
	for(int fromState=1 ; fromState<numFromModelStates ; ++fromState)
	  {
	    compositeHMM.setTransitionProb(fromState+fromModelDelta,0,0);
	  }
      }

  return &compositeHMM;
}



void Application::linkModels(TigrVector<HiddenMarkovModel*> &submodels,
			     int fromModelId,int toModelId,
			     HiddenMarkovModel &compositeHMM,double metaP)
{
  /*
    This method takes every state in the "from" model having a
    transition to zero and redirects that transition so that
    it points to the state(s) that the zero in the "to" model
    can transition to.
   */

  HiddenMarkovModel &fromModel=*submodels[fromModelId];
  HiddenMarkovModel &toModel=*submodels[toModelId];
  int fromModelDelta=submodelOrigins[fromModelId]-1;
  int toModelDelta=submodelOrigins[toModelId]-1;
  int numToModelStates=toModel.countStates();
  int numFromModelStates=fromModel.countStates();

  for(int fromState=1 ; fromState<numFromModelStates ; ++fromState)
    {
      if(!fromModel.doesTransitionExist(fromState,0)) continue;
      double fromP=fromModel.getTransitionProb(fromState,0);
      for(int toState=1 ; toState<numToModelStates ; ++toState)
	{
	  if(!toModel.doesTransitionExist(0,toState)) continue;
	  double toP=toModel.getTransitionProb(0,toState);
	  double newP=fromP*metaP*toP;
	  if(compositeHMM.doesTransitionExist(fromState+fromModelDelta,
				      toState+toModelDelta))
	    throw TigrString("ERROR: transition already exists!");
	  compositeHMM.setTransitionProb(fromState+fromModelDelta,
				    toState+toModelDelta,
				    newP);
	}
      if(toModel.doesTransitionExist(0,0))
	{
	  //cout << "MODEL "<<toModelId<<" HAS 0->0 TRANSITION"<<endl;
	  double newP=fromP*metaP*toModel.getTransitionProb(0,0);
	  closeTransitivity(fromState+fromModelDelta,
			    toModelId,newP,compositeHMM,submodels);
	}
    }
}



void Application::closeTransitivity(int fromState,int intermediateModelId,
				    double fromP,
				    HiddenMarkovModel &compositeHMM,
				    TigrVector<HiddenMarkovModel*> 
				      &submodels)
{
  //cout << "closing transitivity: "<<fromState<<"->"<<intermediateModelId
  //   << "("<<fromP<<")"<<endl;
  for(int toModelId=1 ; toModelId<=numSubmodels ; ++toModelId)
    if(metaModel->doesTransitionExist(intermediateModelId,toModelId))
      {
	double metaP=
	  metaModel->getTransitionProb(intermediateModelId,toModelId);
	HiddenMarkovModel &toModel=*submodels[toModelId];
	int toModelDelta=submodelOrigins[toModelId]-1;
	int numToModelStates=toModel.countStates();
	for(int toState=1 ; toState<numToModelStates ; ++toState)
	  {
	    if(!toModel.doesTransitionExist(0,toState)) continue;
	    double toP=toModel.getTransitionProb(0,toState);
	    double newP=fromP*metaP*toP;
	    if(compositeHMM.doesTransitionExist(fromState,
					toState+toModelDelta))
	      throw TigrString("ERROR: transition already exists!");
	    compositeHMM.setTransitionProb(fromState,
				      toState+toModelDelta,
				      newP);
	  }
	if(toModel.doesTransitionExist(0,0))
	  {
	    //cout << "MODEL "<<toModelId<<" HAS 0->0 TRANSITION"<<endl;
	    double newP=fromP*metaP*toModel.getTransitionProb(0,0);
	    closeTransitivity(fromState,toModelId,newP,compositeHMM,
			      submodels);
	  }
      }
}




