#include <iostream>
#include "TigrMutInfo.H"
using namespace std;



double TigrMutInfo::compute(const TigrContingencyTbl &table)
{
  /*
    Computes mutual information between variables whose joint
    occurrence counts are tabulated in table
   */

  table.computeTotals();

  const int height=table.getSecondDim();
  const int width=table.getFirstDim();
  const double total=(double) table.getGrandTotal();

  TigrJointDistr dist(width,height);
  
  int x, y;
  for(x=0 ; x<width ; ++x)
    for(y=0 ; y<height ; ++y)
      dist[x][y]=table[x][y]/total;

  return compute(dist);
}



double TigrMutInfo::compute(const TigrJointDistr &dist)
{
  /*
    Computes mutual information between two variables whose
    joint distribution is represented in dist
   */

  const int X=table.getFirstDimension();
  const int Y=table.getSecondDimension();
  double mi=0.0;
  
  for(int x=0 ; x<width ; ++x)
    for(int y=0 ; y<height ; ++y)
      {
	double P_sub_xy=dist[x][y];
	double P_sub_x=dist.getMarginalX(x);
	double P_sub_y=dist.getMarginalY(y);
	mi+=P_sub_xy * log(P_sub_xy/(P_sub_x*P_sub_y));
      }

  return mi;
}
