#include <iostream>
#include "TigrChi2FitTest.H"
using namespace std;



TigrChi2FitTest::TigrChi2FitTest(TigrVector<int> &observedCounts,
				     TigrVector<int> &expectedCounts,
				     TigrChi2Table &table)
{
  performTest(observedCounts,expectedCounts,table);
}



bool TigrChi2FitTest::goodFit()
{
  return fit;
}



double TigrChi2FitTest::getChiSquared()
{
  return chiSquared;
}



double TigrChi2FitTest::getP()
{
  return P;
}



void TigrChi2FitTest::performTest(TigrVector<int> &observedCounts,
				    TigrVector<int> &expectedCounts,
				    TigrChi2Table &table)
{
  chiSquared=0;
  int n=observedCounts.size();
  for(int i=0 ; i<n ; ++i)
    {
      int expected=expectedCounts[i];
      if(expected==0) continue;
      int observed=observedCounts[i];
      int numerator=observed-expected;
      chiSquared+=(numerator*numerator/expected);
    }
  int df=n-1;
  P=table.lookupP(df,chiSquared);
  fit=(P>=0.05);
}
