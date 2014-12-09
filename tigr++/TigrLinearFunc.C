#include "TigrLinearFunc.H"
using namespace std;



TigrLinearFunc::TigrLinearFunc(double slope,
			       double intercept,
			       double coefOfDeterm)
  : slope(slope), intercept(intercept),
    coefOfDeterm(coefOfDeterm)
{
   
}



TigrLinearFunc::TigrLinearFunc()
{
   
}



double TigrLinearFunc::operator()(double x) const
{
  return slope*x+intercept;
}



TigrLinearFunc &TigrLinearFunc::operator=(const TigrLinearFunc &f)
{
  slope=f.slope;
  intercept=f.intercept;
  coefOfDeterm=f.coefOfDeterm;
}



double TigrLinearFunc::getCoefDeterm() const
{
  return coefOfDeterm;
}
