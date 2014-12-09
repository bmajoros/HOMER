/**************************************************************
MCE.C : Minimum Classification Error training for pure HMMs
bmajoros@tigr.org 12/9/2004
***************************************************************/
#include <iostream>
#include <fstream>
#include "tigr++/TigrFastaReader.H"
#include "tigr++/Sequence.H"
#include "tigr++/TigrRandom.H"
#include "tigr++/DnaAlphabet.H"
#include "tigr++/TigrCommandLine.H"
#include "tigr++/TigrVector.H"
#include "tigr++/TigrArray1D.H"
#include "HiddenMarkovModel.H"
#include "HigherOrderAlphabet.H"
#include "FastViterbi.H"
using namespace std;


typedef TigrArray2D<int> NGRAM_ARRAY;
typedef TigrArray2D<int>::RowIn2DArray<int> NGRAM_COLUMN;



//======================class HmmNgramCounts============================
struct HmmNgramCounts
{
  // The vector is indexed by state; the map is indexed by ngram ID
  TigrArray1D< TigrMap<int,float> > TP, TN, FP, FN;

  void resize(int);
  void normalize(int numNgrams);
};



//======================class Application============================
class Application
{
  TigrVector<Sequence*> tuningSet1, testSet1;
  TigrVector<Sequence*> tuningSet2, testSet2;
  TigrVector<TigrString> tuningStrings1, testStrings1;
  TigrVector<TigrString> tuningStrings2, testStrings2;
  TigrVector<NGRAM_ARRAY*> tuningNGAs1, testNGAs1;
  TigrVector<NGRAM_ARRAY*> tuningNGAs2, testNGAs2;
  DnaAlphabet alphabet;
  HigherOrderAlphabet *hoa;
  HiddenMarkovModel *hmm1, *hmm2;
  int order, maxIterations, fixedOrder;
  int numStates1, numStates2, alphabetSize;
  FastViterbi *viterbi1, *viterbi2;
  bool useFixedOrder;
  float DEFAULT_UPDATE_COEF;

  void load(const TigrString &filename,TigrVector<Sequence*> &,
	    TigrVector<TigrString> &,TigrVector<NGRAM_ARRAY*> &);
  void MCE();
  void countNgrams(TigrVector<int> &path,
		   NGRAM_ARRAY &ngramArray,
		   TigrArray1D< TigrMap<int,float> > &counts);
  int evaluate(TigrVector<NGRAM_ARRAY*> &posNGAs,
	       FastViterbi &posHmm,
	       FastViterbi &negHmm,
	       HmmNgramCounts &counts1,
	       HmmNgramCounts &counts2);
  float evaluate(TigrVector<NGRAM_ARRAY*> &NGAs1,
		 TigrVector<NGRAM_ARRAY*> &NGAs2,
		 FastViterbi &viterbi1,
		 FastViterbi &viterbi2,
		 HmmNgramCounts &counts1,
		 HmmNgramCounts &counts2);
  void update(HiddenMarkovModel &hmm,
	      HmmNgramCounts &counts,
	      float updateCoef);
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



void HmmNgramCounts::resize(int n)
{
  TP.resize(n);
  TN.resize(n);
  FP.resize(n);
  FN.resize(n);
}



void HmmNgramCounts::normalize(int numNgrams)
{
  int numStates=TP.size();
  for(int state=1 ; state<numStates ; ++state)
    {
      TigrMap<int,float> &TPs=TP[state];
      TigrMap<int,float> &FPs=FP[state];
      TigrMap<int,float> &TNs=TN[state];
      TigrMap<int,float> &FNs=FN[state];

      for(int ngram=0 ; ngram<numNgrams ; ++ngram)
	{
	  float &tp=TPs[ngram];
	  float &tn=TNs[ngram];
	  float &fp=FPs[ngram];
	  float &fn=FNs[ngram];
	  float sum=tp+tn+fp+fn;
	  if(sum>0)
	    {
	      tp/=sum;
	      tn/=sum;
	      fp/=sum;
	      fn/=sum;
	    }
	}
    }
}



void usage()
{
  throw "\n\
MCE <in1.hmm> <in2.hmm> <tuning1.fasta> <tuning2.fasta> <test1.fasta> \n\
    <test2.fasta> <#iterations> <out1.hmm> <out2.hmm> [-s]\n\
\n\
    where: -s <seed> seeds the randomizer,\n\
           -w <weight> specifies initial weight for updates\n\
                       (default: 0.1)\n\
\n";
}



int Application::go(int argc,char *argv[])
{
  // Process command line
  TigrCommandLine cmd(argc,argv,"s:o:w:");
  if(cmd.numArgs()!=9) usage();
  TigrString in1HmmFile=cmd.arg(0);
  TigrString in2HmmFile=cmd.arg(1);
  TigrString tuning1Fasta=cmd.arg(2);
  TigrString tuning2Fasta=cmd.arg(3);
  TigrString test1Fasta=cmd.arg(4);
  TigrString test2Fasta=cmd.arg(5);
  maxIterations=cmd.arg(6).asInt();
  TigrString out1HmmFile=cmd.arg(7);
  TigrString out2HmmFile=cmd.arg(8);
  randomize();
  if(cmd.option('s')) SeedRandomizer(cmd.optParm('s').asUnsigned());
  if(cmd.option('o')) 
    {useFixedOrder=true; fixedOrder=cmd.optParm('o').asInt();}
  else useFixedOrder=false;
  DEFAULT_UPDATE_COEF=0.1;
  if(cmd.option('w'))
    DEFAULT_UPDATE_COEF=cmd.optParm('w').asFloat();

  // Load HMMs
  alphabetSize=alphabet.size();
  hmm1=new HiddenMarkovModel(in1HmmFile,alphabet);
  hmm2=new HiddenMarkovModel(in2HmmFile,alphabet);
  hmm1->useVariableOrder();
  hmm2->useVariableOrder();
  order=hmm1->getOrder();
  if(hmm2->getOrder()!=order) throw "Both HMMs must be of equal order!";
  hoa=new HigherOrderAlphabet(alphabet,order);
  numStates1=hmm1->countStates();
  numStates2=hmm2->countStates();
  viterbi1=new FastViterbi(*hmm1);
  viterbi2=new FastViterbi(*hmm2);

  // Load the training & test sets
  load(tuning1Fasta,tuningSet1,tuningStrings1,tuningNGAs1);
  load(tuning2Fasta,tuningSet2,tuningStrings2,tuningNGAs2);
  load(test1Fasta,testSet1,testStrings1,testNGAs1);
  load(test2Fasta,testSet2,testStrings2,testNGAs2);

  // Perform training
  MCE();

  // Save the results
  hmm1->save(out1HmmFile);
  hmm2->save(out2HmmFile);

  // Clean up
  if(false)
    {
      delete hmm1;
      delete hmm2;
      delete hoa;
      // delete training sequences & NGAs
    }
  return 0;
}



void Application::load(const TigrString &filename,
		       TigrVector<Sequence*> &sequences,
		       TigrVector<TigrString> &strings,
		       TigrVector<NGRAM_ARRAY*> &NGAs)
{
  TigrFastaReader reader(filename);
  TigrString defline, seqStr;
  while(reader.nextSequence(defline,seqStr))
    {
      //### Shorten string
      const int MAX_LEN=25; //30;
      if(seqStr.length()>MAX_LEN) seqStr=seqStr.substr(0,MAX_LEN);

      // Add string to set
      strings.push_back(seqStr);

      // Add sequence to set
      Sequence *seq=new Sequence(seqStr,alphabet);
      sequences.push_back(seq);

      // Add ngramArray to set
      int len=seq->getLength();
      NGRAM_ARRAY *ngramArray=new NGRAM_ARRAY(len,order+1);
      hoa->initNgramArray(*ngramArray,seqStr);
      NGAs.push_back(ngramArray);
    }
}



void Application::countNgrams(TigrVector<int> &path,
			      NGRAM_ARRAY &ngramArray,
			      TigrArray1D< TigrMap<int,float> > &cnts)
{
  int len=path.size();
  for(int pos=0 ; pos<len ; ++pos)
    {
      int state=path[pos];
      if(state<0 || state>=cnts.size()) cout<<state<<" out of bounds"<<endl;
      TigrMap<int,float> &counts=cnts[state];
      int prevNgram=-1;

      // ### Rather than updating all ngrams at this position, we should
      //     probably just update the one actually used by the HMM
      //     (i.e., the longest ngram defined in the emission table)
      NGRAM_COLUMN column=ngramArray[pos];
      for(int o=0 ; o<order ; ++o)
	{
	  int ngram=column[o];
	  if(ngram==prevNgram) break;
	  prevNgram=ngram;
	  if(counts.isDefined(ngram)) counts[ngram]+=1.0;
	  else counts[ngram]=1.0;
	}
    }
}



int Application::evaluate(TigrVector<NGRAM_ARRAY*> &posNGAs,
			  FastViterbi &posHmm,
			  FastViterbi &negHmm,
			  HmmNgramCounts &counts1,
			  HmmNgramCounts &counts2)
{
  double posScore, negScore;
  int n=posNGAs.size(), correct=0;
  for(int i=0 ; i<n ; ++i)
    {
      NGRAM_ARRAY *ngramArray=posNGAs[i];

      TigrVector<int> *posPath=posHmm.getPath(*ngramArray,posScore);
      TigrVector<int> *negPath=negHmm.getPath(*ngramArray,negScore);

      if(posScore>negScore) 
	{
	  countNgrams(*posPath,*ngramArray,counts1.TP);
	  countNgrams(*negPath,*ngramArray,counts2.TN);
	  ++correct;
	}
      else if(posScore<negScore)
	{
	  countNgrams(*posPath,*ngramArray,counts1.FN);
	  countNgrams(*negPath,*ngramArray,counts2.FP);
	}

      delete posPath;
      delete negPath;
    }
  //cout<<"returning correct="<<correct<<endl;//###
  return correct;
}



float Application::evaluate(TigrVector<NGRAM_ARRAY*> &NGAs1,
			    TigrVector<NGRAM_ARRAY*> &NGAs2,
			    FastViterbi &viterbi1,
			    FastViterbi &viterbi2,
			    HmmNgramCounts &counts1,
			    HmmNgramCounts &counts2)
{
  // Initialize the counts arrays
  int numStates1=viterbi1.getHmmGraph().getHMM().countStates();
  int numStates2=viterbi2.getHmmGraph().getHMM().countStates();
  counts1.resize(numStates1);
  counts2.resize(numStates2);

  // Perform classification and take note of misclassified sequences
  int TP=evaluate(NGAs1,viterbi1,viterbi2,counts1,counts2);
  int TN=evaluate(NGAs2,viterbi2,viterbi1,counts2,counts1);
  int FP=NGAs2.size()-TN;
  int FN=NGAs1.size()-TP;
  //cout<<"TP="<<TP<<" TN="<<TN<<" FP="<<FP<<" FN="<<FN<<endl;//###

  // Normalize ngram counts into probabilities
  int numNgrams=hoa->size();
  counts1.normalize(numNgrams);
  counts2.normalize(numNgrams);

  // Return an F-measure of the accuracy
  float sens=TP/float(TP+FN), spec=TP/float(TP+FP);
  float F=2*sens*spec/(sens+spec);
  return F;
}



void Application::update(HiddenMarkovModel &hmm,
			 HmmNgramCounts &counts,
			 float updateCoef)
{
  int numNgrams=hoa->size();
  int numStates=hmm.countStates();
  for(int state=1 ; state<numStates ; ++state)
    {
      TigrMap<int,float> &TPs=counts.TP[state];
      TigrMap<int,float> &TNs=counts.TN[state];
      TigrMap<int,float> &FPs=counts.FP[state];
      TigrMap<int,float> &FNs=counts.FN[state];
      for(int ngram=0 ; ngram<numNgrams ; ++ngram)
	{
	  double oldP=hmm.getEmissionProb(state,ngram);
	  //if(oldP==0.0) continue;
	  //int ngramOrder=hoa->getOrder(ngram);//###
	  //if(ngramOrder<3) continue;//###

	  float tp=TPs[ngram];
	  float tn=TNs[ngram];
	  float fp=FPs[ngram];
	  float fn=FNs[ngram];

	  float update=0;
	  if(fp>tp) update-=3*fp-tp;
	  if(fn>tn) update+=3*fn-tn;
	  
	  //update=0;
	  //if(fp>(tn+tp)) update=-fp-tp;
	  //else if(fn>(tp+tn)) update=fn-tn;


	  if(Random0to1()<0.95) continue;

	  //if(update) cout<<"\t\tupdate="<<update<<endl;
	  if(update<-1) update=-1; else if(update>1) update=1;

	  if(update==0.0) continue;
	  update*=updateCoef;
	  double newP=oldP+update;
	  if(newP<=0) newP=oldP/2;
	  hmm.setEmissionProb(state,ngram,newP);
	}
    }  
}



void Application::MCE()
{
  float currentAccuracy=0;
  for(int iter=1 ; iter<maxIterations ; ++iter)
    {
      // ------------------------------------------------------------------
      // PHASE I -- GETTING A SET OF NGRAMS IMPLICATED IN MISCLASSIFICATION
      // ------------------------------------------------------------------

      // Evaluate classification accuracy & take note of misclassifications
      HmmNgramCounts counts1, counts2;
      float accuracy=evaluate(tuningNGAs1,tuningNGAs2,*viterbi1,*viterbi2,
			      counts1,counts2);
      cout<<iter<<" F="<<accuracy<<" BEST="<<currentAccuracy<<endl;
      currentAccuracy=accuracy;

      // ------------------------------------------------------------------
      // PHASE II -- FINDING THE OPTIMAL WEIGHTING FOR THE UPDATE EQUATIONS
      // ------------------------------------------------------------------
      
      float updateCoef=DEFAULT_UPDATE_COEF;
      for(int rep=0 ; rep<10 ; ++rep)
	{
	  // Make the updates to the HMM parameters
	  HiddenMarkovModel newHmm1(*hmm1);
	  HiddenMarkovModel newHmm2(*hmm2);
	  update(newHmm1,counts1,updateCoef);
	  update(newHmm2,counts2,updateCoef);

	  // Re-evaluate the accuracy
	  HmmNgramCounts counts1, counts2;
	  FastViterbi viterbi1(newHmm1), viterbi2(newHmm2);
	  float accuracy=evaluate(tuningNGAs1,tuningNGAs2,
				  viterbi1,viterbi2,
				  counts1,counts2);
	  cout<<"    NEW F="<<accuracy<<" CURRENT="<<currentAccuracy
	      <<" c="<<updateCoef<<endl;

	  // If accuracy decreased, un-do and try different coefficient
	  if(accuracy>currentAccuracy) 
	    {
	      currentAccuracy=accuracy;
	      *hmm1=newHmm1;
	      *hmm2=newHmm2;
	      break;
	    }
	  else 
	    {
	      if(rep==9)
		{
		  rep=0;
		  updateCoef*=0.5;
		  if(fabs(updateCoef)<0.005) return;
		}
	    }
	}
    }
}



