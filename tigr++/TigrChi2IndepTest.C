using namespace std;
#include <iostream>
#include "TigrChi2IndepTest.H"



TigrChi2IndepTest::TigrChi2IndepTest(TigrArray2D<int> &contingencyTable,
			     TigrChi2Table &distr)
{
   

  runTest(contingencyTable,distr);
}



bool TigrChi2IndepTest::areIndependent()
{
  return independent;
}



double TigrChi2IndepTest::getChiSquared()
{
  return chiSquared;
}



double TigrChi2IndepTest::getP()
{
  return P;
}



void TigrChi2IndepTest::runTest(TigrArray2D<int> &table,TigrChi2Table &distr)
{
  int n1=table.getFirstDim(), n2=table.getSecondDim();
  TigrIntArray1D rowSums(n1), columnSums(n2);
  rowSums.setAllTo(0);
  columnSums.setAllTo(0);
  float tableSum=0;

  for(int i=0 ; i<n1 ; ++i)
    for(int j=0 ; j<n2 ; ++j)
      {
	int x=table[i][j];
	rowSums[i]+=x;
	columnSums[j]+=x;
	tableSum+=x;
      }

  chiSquared=0;
  for(int i=0 ; i<n1 ; ++i)
    for(int j=0 ; j<n2 ; ++j)
      {
	int e=int(rowSums[i]/tableSum*columnSums[j]);
	if(e==0) continue;
	int o=table[i][j];
	int diff=o-e;
	float term=diff*diff/float(e);
	 
	chiSquared+=term;
      }
  
  int df=(n1-1)*(n2-1);
  P=distr.lookupP(df,chiSquared);

  independent=(P>0.05);
}
