/**************************************************************
hmm-substitute : for substituting the states of one HMM with
                 those of another
bmajoros@tigr.org 12/11/2004
***************************************************************/
#include <iostream>
#include <fstream>
#include "tigr++/DnaAlphabet.H"
#include "tigr++/TigrCommandLine.H"
#include "HiddenMarkovModel.H"
#include "HigherOrderAlphabet.H"
using namespace std;


//======================class Application============================
class Application
{
  DnaAlphabet alphabet;
  HigherOrderAlphabet *hoa;
  HiddenMarkovModel *bigHmm, *littleHmm;
  int order;
  int numBigStates, numLittleStates, alphabetSize;

  void sub(int stateInBigHmm);
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
hmm-substitute <big.hmm> <little.hmm> <state-in-big-hmm> <out.hmm>\n\
\n";
}



int Application::go(int argc,char *argv[])
{
  // Process command line
  TigrCommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=4) usage();
  TigrString bigHmmFile=cmd.arg(0);
  TigrString littleHmmFile=cmd.arg(1);
  int stateInBigHmm=cmd.arg(2).asInt();
  TigrString outHmmFile=cmd.arg(3);

  // Load HMMs
  alphabetSize=alphabet.size();
  bigHmm=new HiddenMarkovModel(bigHmmFile,alphabet);
  littleHmm=new HiddenMarkovModel(littleHmmFile,alphabet);
  order=bigHmm->getOrder();
  if(littleHmm->getOrder()!=order) 
    throw "Both HMMs must be of equal order!";
  hoa=new HigherOrderAlphabet(alphabet,order);
  numBigStates=bigHmm->countStates();
  numLittleStates=littleHmm->countStates();

  // Perform substitution
  sub(stateInBigHmm);

  // Save the result
  bigHmm->save(outHmmFile);

  // Clean up
  if(false)
    {
      delete bigHmm;
      delete littleHmm;
      delete hoa;
    }
  return 0;
}



void Application::sub(int startState)
{
  int numNgrams=hoa->size();
  int bigState=startState;
  for(int littleState=1 ; littleState<numLittleStates ; ++littleState)
    {
      cout<<littleState<<" -> "<<bigState<<endl;
      for(int ngram=0 ; ngram<numNgrams ; ++ngram)
	{
	  double p=littleHmm->getEmissionProb(littleState,ngram);
	  bigHmm->setEmissionProb(bigState,ngram,p);
	}
      ++bigState;
    }
}



