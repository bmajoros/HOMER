#include <iostream>
#include "TigrPipe.H"
using namespace std;



TigrPipe::TigrPipe(const TigrString &command,const TigrString &mode)
  : TigrFile()
{
  fp=popen(command.c_str(),mode.c_str());
  if(!fp) 
    throw TigrString("Can't open pipe for command \"")+
      command+"\"";
  filename=command;
  this->mode=mode;
}



void TigrPipe::close()
{
  if(fp!=NULL)
    pclose(fp);
  fp=NULL;
}
