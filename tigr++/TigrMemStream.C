#include <iostream>
#include "TigrMemStream.H"
using namespace std;



TigrMemStream::TigrMemStream(const char *p)
  : p(p)
{
   
}



TigrString TigrMemStream::readLine()
{
  const char *q=p;
  while(true)
    switch(*q)
      {
      case '\n': 
      case '\0': 
	goto end;
      default:
	++q;
      }
 end:
  const char *oldP=p;
  unsigned length=q-p;
  p=q+1;
  return TigrString(oldP,length);
}



bool TigrMemStream::eof()
{
  return *p=='\0';
}
