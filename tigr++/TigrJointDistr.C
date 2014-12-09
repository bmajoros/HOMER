#include <iostream>
#include "TigrJointDistr.H"
using namespace std;



TigrJointDistr::TigrJointDistr(int X,int Y)
  : TigrArray2D<double>(X,Y), marginalX(X), marginalY(Y)
{
   
}



double TigrJointDistr::getMarginalX(int x)
{
  return marginalX[x];
}



double TigrJointDistr::getMarginalY(int y)
{
  return marginalY[y];
}



void TigrJointDistr::computeMarginals()
{
  marginalX.setAllTo(0.0);
  marginalY.setAllTo(0.0);

  const int X=getFirstDimension();
  const int Y=getSecondDimension();
  int x, y;
  for(x=0 ; x<X ; ++x)
    for(y=0 ; y<Y ; ++y)
      {
	double entry=(*this)[x][y];
	marginalX[x]+=entry;
	marginalY[Y]+=entry;
      }  
}
