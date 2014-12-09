#include "Factorial.C"
#include <iostream>

double logFactorial(int x);


int main(void)
{
  for(int i=3 ; i<100 ; ++i)
    {
      cout << logFactorial(i) << " vs. " << Factorial::global.logFactorial(i) << endl;
    }
  return 0;
}


double logFactorial(int x)
{
  double sum=0;
  for(int i=2 ; i<=x ; ++i)
    sum+=log((long double)i);
  return sum;
}
