#include <iostream>
#include "NormalDistribution.C"
#include "Factorial.C"

int main(void)
{
  NormalDistribution p(80,10);
  for(int i=0 ; i<=100 ; ++i)
    cout<<i<<" "<<p.density(i)<<endl;
  return 0;
}

