#include <math.h>
#include <iostream>
#include "TigrSummaryStats.H"
using namespace std;


TigrSummaryStats::TigrSummaryStats(const TigrVector<int> &v)
{
   

  compute(v);
}



TigrSummaryStats::TigrSummaryStats(const TigrVector<float> &v)
{
   

  compute(v);
}



TigrSummaryStats::TigrSummaryStats(const TigrVector<double> &v)
{
   

  compute(v);
}



double TigrSummaryStats::getMax()
{
  return max;
}



double TigrSummaryStats::getMean()
{
  return mean;
}



double TigrSummaryStats::getMin()
{
  return min;
}



double TigrSummaryStats::getStdDev()
{
  return stddev;
}



int TigrSummaryStats::getN()
{
  return n;
}



void TigrSummaryStats::compute(const TigrVector<double> &v)
{
  n=v.size();
  double sumX=0.0, sumXX=0.0;
  int i;
  for(i=0 ; i<n ; ++i)
    {
      double x=v[i];
      sumX+=x;
      sumXX+=x*x;
      if(i==0) min=max=x;
      else if(x<min) min=x;
      else if(x>max) max=x;
    }
  mean=sumX/n;
  stddev=sqrt((sumXX-sumX*sumX/n)/(n-1.0));
}



void TigrSummaryStats::compute(const TigrVector<float> &v)
{
  n=v.size();
  double sumX=0.0, sumXX=0.0;
  int i;
  for(i=0 ; i<n ; ++i)
    {
      double x=(double) v[i];
      sumX+=x;
      sumXX+=x*x;
      if(i==0) min=max=x;
      else if(x<min) min=x;
      else if(x>max) max=x;
    }
  mean=sumX/n;
  stddev=sqrt((sumXX-sumX*sumX/n)/(n-1.0));
}



void TigrSummaryStats::compute(const TigrVector<int> &v)
{
  n=v.size();
  double sumX=0.0, sumXX=0.0;
  int i;
  for(i=0 ; i<n ; ++i)
    {
      double x=(double) v[i];
      sumX+=x;
      sumXX+=x*x;
      if(i==0) min=max=x;
      else if(x<min) min=x;
      else if(x>max) max=x;
    }
  mean=sumX/n;
  stddev=sqrt((sumXX-sumX*sumX/n)/(n-1.0));
}
