#include <iostream>
#include "TigrFastaReader.H"
using namespace std;



TigrRegex TigrFastaReader::deflineRegex("^\\s*>\\s*(\\S+)(.*)");



TigrFastaReader::TigrFastaReader(const TigrString &filename,
				 Alphabet &alphabet)
  : file(filename), alphabet(alphabet)
{
}



TigrVector< pair<TigrString,TigrString> > *TigrFastaReader::readAll()
{
  TigrVector< pair<TigrString,TigrString> > *v=new 
    TigrVector< pair<TigrString,TigrString> >;
  TigrString defline, sequence;
  while(nextSequence(defline,sequence))
    v->push_back(pair<TigrString,TigrString>(defline,sequence));
  return v;
}



bool TigrFastaReader::nextSequence(TigrString &defline,TigrString &sequence)
{
   
  if(file.eof()) return false;
  if(cache.length()>0) 
    {
      defline=cache;
      cache="";
    }
  else defline=file.readLine();
  if(file.eof()) return false;

   
  sequence="";
  while(true)
    {
      TigrString line=file.readLine();
      if(file.eof()) break;
      if(line[0]=='>')
	{
	  cache=line;
	  break;
	}
      line.trimWhitespace();
      sequence+=line;
    }
  sequence.toupper();
  maskStrangeChars(sequence);
  return true;
}



void TigrFastaReader::maskStrangeChars(TigrString &s)
{
  int len=s.length();
  for(int i=0 ; i<len ; ++i)
    if(!alphabet.isDefined(s[i])) s[i]='N';
}



void TigrFastaReader::parseDefline(const TigrString &defline,TigrString &id,
				   TigrString &remainder)
{
  if(!deflineRegex.match(defline)) 
    throw TigrString("TigrFastaReader::parseDefline() failed to parse: ")+
      defline;
  id=deflineRegex[1];
  remainder=deflineRegex[2];
}



TigrString TigrFastaReader::getId(const TigrString &defline)
{
  TigrString id, remainder;
  parseDefline(defline,id,remainder);
  return id;
}



