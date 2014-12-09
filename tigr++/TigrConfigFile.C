#include <stdlib.h>
#include <fstream>
#include "TigrFile.H"
#include "TigrConfigFile.H"
using namespace std;



TigrConfigFile::TigrConfigFile(TigrString filename)
{
  load(filename);
}



TigrConfigFile::TigrConfigFile()
{
}



TigrString TigrConfigFile::lookup(TigrString attr)
{
  return dict[attr];
}



TigrString TigrConfigFile::lookupOrDie(TigrString attr)
{
  return charLookupOrDie(attr);
}



const char *TigrConfigFile::charLookupOrDie(TigrString attr)
{
  if(dict.find(attr)==dict.end()) 
    throw attr+" is not defined in the configuration file";
  return dict[attr].c_str();
}



double TigrConfigFile::getDoubleOrDie(TigrString attr)
{
  return atof(charLookupOrDie(attr));
}



float TigrConfigFile::getFloatOrDie(TigrString attr)
{
  return (float) atof(charLookupOrDie(attr));
}



int TigrConfigFile::getIntOrDie(TigrString attr)
{
  return atoi(charLookupOrDie(attr));
}



bool TigrConfigFile::getBoolOrDie(TigrString attr)
{
  TigrString b=charLookupOrDie(attr);
  if(b=="true" || b=="t" || b=="yes" || b=="y") return true;
  if(b=="false" || b=="f" || b=="no" || b=="n") return false;
  throw attr+" must be true or false in config file";
}



long TigrConfigFile::getLongOrDie(TigrString attr)
{
  return atol(charLookupOrDie(attr));
}



void TigrConfigFile::enter(TigrString key,TigrString value)
{
  dict[key]=value;
}



void TigrConfigFile::load(TigrString fname)
{
  TigrFile file(fname);
  if(!file.isOpen())
    throw TigrString("Error: Can't open file ")+fname+
      " in TigrConfigFile::load()";
  while(!file.eof())
    {
      TigrString nextline=file.getline();
      
      if(file.eof()) break;
      const char *p=nextline.c_str();
      TigrStrTokenizer parser(p,"=\n\t \r");
      processLine(parser);
    }
}



void TigrConfigFile::processLine(TigrStrTokenizer &parser)
{
  if(!parser.hasMoreTokens()) return;  
  TigrString key=parser.nextToken();
  if('#'==key[0]) return; // comment
  if(!parser.hasMoreTokens())  
    throw TigrString("Syntax error in configuration file: "+key);
  const char *value=parser.nextToken();
  dict[key]=value; 
}
