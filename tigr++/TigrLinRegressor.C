#include <math>
#include "TigrString.H"
#include "TigrStacktrace.H"
#include "TigrLinRegressor.H"
using namespace std;



TigrLinearFunc TigrLinRegressor::compute(double sumX,double sumXX,
					 double sumY,double sumYY,
					 double sumXY,int n)
{
  const double yBar=sumY/n;
  const double xBar=sumX/n;

  const double Sxx=sumXX-sumX*sumX/n;
  if(Sxx==0.0) throw TigrString("Error in TigrLinRegressor: Sxx is 0");

  const double Syy=sumYY-sumY*sumY/n;
  const double Sxy=sumXY-sumX*sumY/n;

  const double slope=Sxy/Sxx;
  const double intercept=yBar-xBar*slope;
   
  const double r=Sxy/sqrt(Syy*Sxx);
  const double coefDetermination=r*r;

  return TigrLinearFunc(slope,intercept,coefDetermination);
}



TigrLinearFunc TigrLinRegressor::regress(DblPointVector &points)
{
  int numPoints=points.size();
  if(numPoints<2) 
    throw TigrString("TigrLinRegressor::regress(): "
		     "Regression requires at least 2 points");

  double sumXX=0.0, sumYY=0.0, sumXY=0.0, sumY=0.0, sumX=0.0;

  for(int i=0 ; i<numPoints ; ++i)
    {
      DblPoint &point=points[i];
      double x=point.first;
      double y=point.second;

      sumX+=x;
      sumY+=y;
      sumXX+=x*x;
      sumYY+=y*y;
      sumXY+=x*y;
    }

  return compute(sumX,sumXX,sumY,sumYY,sumXY,numPoints);
}



TigrLinearFunc TigrLinRegressor::regress(TigrVector<double> &yValues)
{
  int numPoints=yValues.size();
  if(numPoints<2) 
    throw TigrString("Regressing on <2 points");

  double sumX=0.0, sumY=0.0, sumXX=0.0, sumYY=0.0, sumXY=0.0;
  for(int i=0 ; i<numPoints ; ++i)
    {
      double x=(double) i;
      double y=yValues[i];

      sumX+=x;
      sumY+=y;
      sumXX+=x*x;
      sumYY+=y*y;
      sumXY+=x*y;
    }

  return compute(sumX,sumXX,sumY,sumYY,sumXY,numPoints);
}
