#include <iostream>
#include "ExtremeValueDistribution.C"
#include "Factorial.C"

int main(void)
{
  ExtremeValueDistribution p(47.0,2.0);
  for(int i=0 ; i<=100 ; ++i)
    cout<<i<<" "<<p.density(i)<<endl;
  return 0;
}

