
// viterbi-trainer.C

#include <iostream>
#include "tigr++/TigrCommandLine.H"
#include "tigr++/TigrVector.H"
#include "tigr++/Sequence.H"
#include "tigr++/DnaAlphabet.H"
#include "tigr++/TigrArray2D.H"
#include "tigr++/TigrFastaReader.H"
#include "HigherOrderAlphabet.H"
#include "FastViterbi.H"
using namespace std;


class Application
{
  TigrVector<Sequence*> trainingSet;
  TigrVector<TigrString> trainingSequences;
  DnaAlphabet alphabet;
  HigherOrderAlphabet *oldHOA, *newHOA;

  void load(TigrString filename);
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
    if(cmd.numArgs()!=3)
      throw string("viterbi-trainer <template.hmm> <training-features.fasta> <order>");
    TigrString inHmmFile=cmd.arg(0);
    TigrString inTrainFile=cmd.arg(1);
    int newOrder=cmd.arg(2).asInt();
    //TigrString outHmmFile=cmd.arg(3);

    // Load HMM and perform some initialization
    HiddenMarkovModel oldHmm(inHmmFile,alphabet);
    FastViterbi viterbi(oldHmm);
    int oldOrder=oldHmm.getOrder();
    oldHOA=new HigherOrderAlphabet(alphabet,oldOrder);

    // Load training set
    load(inTrainFile);
    int numTrain=trainingSet.size();
    newHOA=new HigherOrderAlphabet(alphabet,newOrder);

    // Iterate through all training sequences
    for(int i=0 ; i<numTrain ; ++i)
      {
	Sequence &seq=*trainingSet[i];
	const TigrString &seqStr=trainingSequences[i];
	TigrArray2D<int> ngramArray(seqStr.length(),oldOrder+1);
	oldHOA->initNgramArray(ngramArray,seqStr);
	double pathScore;
	vector<int> &path=*viterbi.getPath(seq,ngramArray,pathScore);
	int len=path.size();
	for(int j=0 ; j<len ; ++j)
	  {
	    int begin=j-newOrder;
	    if(begin<0) begin=0;
	    int ngramLen=j-begin+1;
	    TigrString ngram=seqStr.substr(begin,ngramLen);
	    cout<<ngram<<" "<<path[j]<<endl;
	  }
      }
    
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

