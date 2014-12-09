#include <ctype.h>
#include "TigrString.H"
#include "TigrStrTokenizer.H"
using namespace std;


TigrStrTokenizer::TigrStrTokenizer(const char *wspace)
  : p(NULL), source(NULL), wspace(wspace)
{
}



TigrStrTokenizer::TigrStrTokenizer(const char *src,
				   const char *wspace)
  : wspace(wspace), source(src), p(src)
{
}



TigrStrTokenizer::TigrStrTokenizer(const TigrString &src,
				   const char *wspace)
  : wspace(wspace), source(src.c_str()), p(src.c_str())
{
}



TigrVector<TigrString*> *TigrStrTokenizer::getAllTokens()
{
  TigrVector<TigrString*> &toks=
    *new TigrVector<TigrString*>;
  while(hasMoreTokens())
    {
      TigrString *s=new TigrString(nextToken());
      toks.push_back(s);
    }
  return &toks;
}



TigrVector<TigrString*> *TigrStrTokenizer::tokenize(const char *src,
						    const char *wspace)
{
  TigrStrTokenizer tokenizer(src,wspace);
  TigrVector<TigrString*> *toks=tokenizer.getAllTokens();
  return toks;
}



TigrVector<TigrString> *TigrStrTokenizer::getTokenStrings()
{
  TigrVector<TigrString> &toks=
    *new TigrVector<TigrString>;
  while(hasMoreTokens())
    toks.push_back(TigrString(nextToken()));
  return &toks;
}



bool TigrStrTokenizer::hasMoreTokens()
{
  skipWhiteSpace();
  return p && (*p!='\0');
}



bool TigrStrTokenizer::isWhiteSpace(char c)
{
  if(!*wspace) return !isalpha(c);

  const char *wptr=wspace;
  while(*wptr)
    {
      if(*wptr==c) return true;
      ++wptr;
    }
  return false;
}



const char *TigrStrTokenizer::nextToken()
{
  skipWhiteSpace();
  buffer="";
  while(*p && !isWhiteSpace(*p))
    {
      buffer+=*p;
      ++p;
    }
  const char *rval=buffer.c_str();
  return rval;
}



StrVectPair* TigrStrTokenizer::getTokensAndSeparators(
  const TigrString &src)
{
  WhereType state=AT_START;
  TigrString spc, tok;
  StrVectPair &thePair=*new StrVectPair;
  for(const char *ptr=src.c_str() ; *ptr ; ++ptr)
    {
      if(!isWhiteSpace(*ptr))
	{
	  if(state==IN_SPACE)
	    {
	      thePair.second.push_back(spc);
	      spc="";
	    }
	  tok+=*ptr;
	  state=IN_TOKEN;
	}
      else
	{
	  if(state==IN_TOKEN)
	    {
	      thePair.first.push_back(tok);
	      tok="";
	    }
	  spc+=*ptr;
	  state=IN_SPACE;
	}
    }

  if(tok.length()) thePair.first.push_back(tok);
  if(spc.length()) thePair.second.push_back(spc);

  return &thePair;
}



void TigrStrTokenizer::skipWhiteSpace()
{
  while(*p && isWhiteSpace(*p)) p++;
}



void TigrStrTokenizer::tokenizeThis(const TigrString &str)
{
  p=str.c_str();
  source=str.c_str();
}



void TigrStrTokenizer::tokenizeThis(const char *str)
{
  p=str;
  source=str;
}
