/**************************************************************
baum-welch.C : For training an HMM using Expectation Maximization
bmajoros@tigr.org 1/1/2003

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/

#include <string>
#include <iostream>
#include <fstream>
#include "tigr++/TigrFastaReader.H"
#include "tigr++/Sequence.H"
#include "tigr++/TigrRandom.H"
#include "tigr++/DnaAlphabet.H"
#include "tigr++/TigrCommandLine.H"
#include "tigr++/TigrVector.H"
#include "HMMreader.H"
#include "HMMbuilder.H"
#include "BaumWelch.H"
#include "HigherOrderAlphabet.H"
using namespace std;


class Application
{
  TigrVector<Sequence*> trainingSet;
  TigrVector<TigrString> trainingSequences;
  DnaAlphabet alphabet;

  void load(TigrString);
public:
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
  catch(...)
    {
      cerr << "unknown exception caught in main()" << endl;
      return -1;
    }
  return 0;
}



void usage()
{
  throw "\n\
baum-welch <structure.hmms> <examples.fasta> <#iterations> \n\
           <outfile> <min-sample-size> <order> [-v -s]\n\
where -v = use variable order throughout (numerically unstable!)\n\
      -s <seed> = seed the random number generator\n\
";
}



int Application::go(int argc,char *argv[])
{
  // Process command line
  TigrCommandLine cmd(argc,argv,"nvs:");
  if(cmd.numArgs()!=6) usage();
  TigrString structureFile=cmd.arg(0);
  TigrString trainFilename=cmd.arg(1);
  int maxIterations=cmd.arg(2).asInt();
  TigrString outfile=cmd.arg(3);
  int minSampleSize=cmd.arg(4).asInt();
  int order=cmd.arg(5).asInt();
  bool variableOrder=cmd.option('v');
  randomize();
  if(cmd.option('s')) SeedRandomizer(cmd.optParm('s').asUnsigned());
  bool useNormalized=cmd.option('n');

  // Load HMM structure
  HigherOrderAlphabet higherOrderAlphabet(alphabet,order);
  HMMreader reader(alphabet,higherOrderAlphabet);
  HiddenMarkovModel *hmm=reader.read(structureFile,order);
  if(variableOrder) hmm->useVariableOrder();

  // Load the training set
  load(trainFilename);

  // Train HMM using the "Higher Order Normalized Baum-Welch" algorithm
  BaumWelch(*hmm,maxIterations,trainingSet,trainingSequences,
	    minSampleSize,order,higherOrderAlphabet,useNormalized);

  // Save the results
  hmm->save(outfile);

  return 0;
}



void Application::load(TigrString filename)
{
  TigrFastaReader reader(filename);
  TigrString defline, sequence;
  while(reader.nextSequence(defline,sequence))
    {
      trainingSequences.push_back(sequence);
      Sequence *seq=new Sequence(sequence,alphabet);
      trainingSet.push_back(seq);
    }
}






