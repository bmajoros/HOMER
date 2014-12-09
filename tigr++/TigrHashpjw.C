#include <iostream>
#include "TigrHashpjw.H"
using namespace std;



unsigned TigrHashpjw(const char *s,unsigned tableSize)
{
  return TigrHashpjw(TigrString(s),tableSize);
}



unsigned TigrHashpjw(const TigrString &s,unsigned tableSize)
{
  return TigrHashpjw(TigrString(s),tableSize);
}



unsigned TigrHashpjw(TigrString s,unsigned tableSize)
{
  unsigned length=s.length();
  if(length<5) 
    {
      s=s+s;
      length*=2;
    }
  int h=0;
  for(int i=0 ; i<length ; ++i)
    {
      h=(h<<4)+s[i];
      int g=h & 0xf000;
      if(g) h=h^(g>>8);
    }
  return h % tableSize;
}
