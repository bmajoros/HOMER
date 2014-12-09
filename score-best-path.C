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


class Application
{
  DnaAlphabet alphabet;
  HigherOrderAlphabet *higherOrderAlphabet;
  int numSymbols;
  int order;

  HiddenMarkovModel *loadHMM(const TigrString &);
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
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    TigrCommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=2)
      throw string("score-best-path <*.hmm> <*.fasta>\n");
    TigrString hmmFilename=cmd.arg(0);
    TigrString fastaFilename=cmd.arg(1);

    HiddenMarkovModel &hmm=*new HiddenMarkovModel(hmmFilename,alphabet);
    FastViterbi viterbi(hmm);
    order=hmm.getOrder();
    higherOrderAlphabet=new HigherOrderAlphabet(alphabet,order);

    TigrFastaReader fastaReader(fastaFilename);
    TigrString defline, sequence;
    while(fastaReader.nextSequence(defline,sequence))
      {
	Sequence seq(sequence,alphabet);
	TigrArray2D<int> ngramArray(seq.getLength(),order+1);
	higherOrderAlphabet->initNgramArray(ngramArray,sequence);
	double pathScore;
	vector<int> &path=*viterbi.getPath(seq,ngramArray,pathScore);
	cout<<pathScore<<endl;
      }
    return 0;
  }



HiddenMarkovModel *Application::loadHMM(const TigrString &filename)
{
  HiddenMarkovModel *H=new HiddenMarkovModel(filename,alphabet);
  return H;
}



