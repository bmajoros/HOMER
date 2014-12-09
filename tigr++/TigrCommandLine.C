#include <unistd.h>
#include "TigrCommandLine.H"
using namespace std;



TigrCommandLine::TigrCommandLine(int argc,char *argv[],const char *options)
  : argc(argc), argv(argv)
{
  loopThroughOptions(options);
}



TigrString TigrCommandLine::arg(int i)
{
  return args[i];
}



TigrString TigrCommandLine::optParm(char c)
{
  if(optionParm.find(c)==optionParm.end()) 
    return TigrString("");
  return optionParm[c];
}



bool TigrCommandLine::option(char c)
{
  return usedOption.find(c)==usedOption.end() ? false : usedOption[c];
}



bool TigrCommandLine::takesParameter(char c,const char *options)
{
  for(int i=0 ; options[i] ; ++i)
    if(options[i]==c)
      return options[i+1]==':';
  return false;
}



int TigrCommandLine::numArgs()
{
  return args.size();
}



void TigrCommandLine::loopThroughOptions(const char *options)
{
  int c;
  while((c=getopt(argc,argv,options))!=EOF)
    {
      usedOption[c]=true;
      if(takesParameter(c,options))
	optionParm[c]=TigrString(optarg);
    }

  for(int i=optind ; i<argc ; ++i)
    args.push_back(argv[i]);
}
