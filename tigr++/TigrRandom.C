#include <time.h>
#include <stdlib.h>
#include <limits.h>
#include <iostream>
#include "TigrRandom.H"
using namespace std;



float Random0to1()
{
  const int randInt=RandomNumber(LARGEST_RANDOM_NUMBER);

  return randInt/(float)(LARGEST_RANDOM_NUMBER-1);
}



float RandomFloat(float from,float to)
{
  return RandomFloat(to-from)+from;
}



float RandomFloat(float range)
{
  return Random0to1()*range;
}



float RandomGaussian(float min,float max,int n)
{
  float sum=0.0;
  int i;
  for(i=0 ; i<n ; ++i)
    sum+=RandomFloat(min,max);
  return sum/n;
}



int RandomNumber(int n)
{
  return n>1 ? rand() % n : 0;
}



unsigned GetRandomSeed()
{
  randomize();  
  return unsigned(RandomNumber(INT_MAX));
}



void SeedRandomizer(unsigned s)
{
  srand(s);
}



void randomize()
{
  unsigned seed=(unsigned) time(0);
  srand(seed);
}
