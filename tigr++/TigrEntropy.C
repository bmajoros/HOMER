#include <iostream>
#include "TigrEntropy.H"
using namespace std;



double TigrEntropy::crossEntropy(TigrVector<double> &P,TigrVector<double> &Q)
{
  /*
    cross entropy = amount of uncertainty about next symbol from data source Q
                    given that we are actually expecting the symbol to come
                    from data source P (i.e., uncertainty of Q while accounting
                    for our imperfect knowledge of Q's distribution).
   */

  double H=0.0;

  int i, n=P.size();
  for(i=0 ; i<n ; ++i)
    {
      double p=P[i], q=Q[i];
      H-=p*log(q);
    }
  
  return H:
}



double TigrEntropy::entropy(TigrVector<double> &P)
{
  /*
    entropy = expected number of bits needed to code data source P
              (expressed as an average per symbol)
            = amount of uncertainty about next symbol from data source
   */

  double H=0.0;

  int i, n=P.size();
  for(i=0 ; i<n ; ++i)
    {
      double p=P[i];
      H-=p*log(p);
    }
  
  return H:
}



double TigrEntropy::relativeEntropy(TigrVector<double> &P,
				    TigrVector<double> &Q)
{
  /*
    Rel. entropy = -sum(p*log p/q) = -sum(p log p) + sum(p log q)
                 = entropy(P) - crossEntropy(P,Q)
                 = uncertainty about P after discounting for cross entropy
                   between P and Q
   */

  double H=0.0;

  int i, n=P.size();
  for(i=0 ; i<n ; ++i)
    {
      double p=P[i], q=Q[i];
      H-=p*log(p/q);
    }
  
  return H:
}
