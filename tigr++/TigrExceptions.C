#include <iostream>
#include "TigrExceptions.H"
using namespace std;



ArrayIndexException::ArrayIndexException(long index,const TigrString &msg)
  : RootException("")
{
   

  installMessage(index,msg);
}



RootException::RootException(const TigrString &s)
  : message(s)
{
   
}



RootException::RootException(const char *cp)
  : message(cp)
{
   
}



const TigrString &RootException::getMessage() const
{
  return message;
}



void ArrayIndexException::installMessage(long index,const TigrString &msg)
{
  message=TigrString("Invalid index (")+index+") "+msg;
}
