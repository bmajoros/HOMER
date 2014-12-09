#include <iostream>
#include "PoissonDistribution.C"
#include "Factorial.C"

int main(void)
{
  PoissonDistribution p(25);
  for(int i=0 ; i<=100 ; ++i)
    cout<<i<<" "<<p.density(i)<<endl;
  return 0;
}

