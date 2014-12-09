#include <iostream>
#include "TigrRegex.H"
using namespace std;



TigrRegex::TigrRegex(const TigrString &regex)
  : regex(preprocess(regex))
{
  reg.start=NULL;
  reg.end=NULL;

  compile();
}



TigrRegex::~TigrRegex()
{
  regfree(&patbuf);
  free(reg.start);
  free(reg.end);
}



TigrString TigrRegex::operator[](int rNum)
{
  int fr=reg.start[rNum], to=reg.end[rNum];
  int length=to-fr;
  if(!substrate) throw "No substrate in TigrRegex::operator[]";
  return substrate->substring(fr,length);
}



TigrString TigrRegex::getEntireMatch() 
{
  return (*this)[0];
}



TigrString TigrRegex::preprocess(const TigrString &express)
{
  TigrString rval;
  const char *ptr=express.c_str();
  while(*ptr)
    {
      if(*ptr=='\\')
	{
	  ++ptr;
	  switch(*ptr)
	    {
	    case 'P':  rval+="[^[:print:]]"; break;
	    case 'U':  rval+="[^[:upper:]]"; break;
	    case 'a':  rval+="[[:alpha:]]"; break;
	    case 's':  rval+="[[:blank:][:space:]]"; break;
	    case 'S':  rval+="[^[:blank:][:space:]]"; break;
	    case 'G':  rval+="[^[:graph:]]"; break;
	    case 'L':  rval+="[^[:lower:]]"; break;
	    case 'u':  rval+="[[:upper:]]"; break;
	    case 'c':  rval+="[[:cntrl:]]"; break;
	    case 'D':  rval+="[^[:digit:]]"; break;
	    case 'x':  rval+="[[:xdigit:]]"; break;
	    case 'w':  rval+="[[:alnum:]]"; break;
	    case 'l':  rval+="[[:lower:]]"; break;
	    case 'C':  rval+="[^[:cntrl:]]"; break;
	    case '\0': rval+='\\'; continue;
	    case 'd':  rval+="[[:digit:]]"; break;
	    case 'X':  rval+="[^[:xdigit:]]"; break;
	    case 'p':  rval+="[[:print:]]"; break;
	    case 'W':  rval+="[^[:alnum:]]"; break;
	    case 'A':  rval+="[^[:alpha:]]"; break;
	    case 'g':  rval+="[[:graph:]]"; break;
	    }
	  ++ptr;
	}
      else
	{
	  rval+=*ptr;
	  ++ptr;
	}
    }
  return rval;
}



TigrString TigrRegex::substitute(const TigrString &regex,
				 const TigrString &replacement,
				 const TigrString &substrate)
{
  TigrRegex r(regex);
  return r.substitute(substrate,replacement);
}



TigrString TigrRegex::substitute(const TigrString &s,
				 const TigrString &replacement)
{
  TigrString rval;
  int len=s.length(), here=0;
  while(here<len)
    {
      int matchIndex=re_search(&patbuf,s.c_str(),s.size(),here,s.size(),&reg);
      if(matchIndex<0) break;
      rval+=s.substring(here,matchIndex-here);
      rval+=replacement;
      int matchLen=reg.end[0]-reg.start[0];
      if(matchLen==0) ++matchLen;
      here=matchIndex+matchLen;
    }
  if(here<len)
    rval+=s.substring(here,len-here);
  return rval;
}



TigrVector<TigrString*> *TigrRegex::split(const TigrString &regex,
					  const TigrString &substrate)
{
  TigrRegex r(regex);
  return r.split(substrate);
}



TigrVector<TigrString*> *TigrRegex::split(const TigrString &substrate)
{
  TigrVector<TigrString*> *fields=new TigrVector<TigrString*>;

  int len=substrate.length(), here=0;
  while(here<len)
    {
      int nextMatchIndex=
	re_search(&patbuf,substrate.c_str(),substrate.size(),here,
		  substrate.size(),&reg);
      if(nextMatchIndex<0) break;
      if(nextMatchIndex>here)
	fields->push_back(
	  new TigrString(substrate.substring(here,nextMatchIndex-here)));
      int matchLen=
	reg.end[0]-reg.start[0];
      if(matchLen==0) 
	++matchLen;
      here=nextMatchIndex+matchLen;
    }
  if(here<len)
    fields->push_back(
      new TigrString(substrate.substring(here,len-here)));

  return fields;
}



bool TigrRegex::match(const TigrString &regex,const TigrString &substrate)
{
  TigrRegex r(regex);
  return r.match(substrate);
}



bool TigrRegex::match(const TigrString &s)
{
  int matchLength=re_match(&patbuf,s.c_str(),s.size(),0,&reg);
  if(matchLength<0)
    {
      substrate=NULL;
      return false;
    }
  substrate=&s;
  return true;
}



bool TigrRegex::search(const TigrString &regex,const TigrString &substrate)
{
  TigrRegex r(regex);
  return r.search(substrate);
}



bool TigrRegex::search(const TigrString &s)
{
  int matchIndex=re_search(&patbuf,s.c_str(),s.size(),0,s.size(),&reg);
  if(matchIndex<0)
    {
      substrate=NULL;
      return false;
    }
  substrate=&s;
  return true;
}



int TigrRegex::getNumSubexpressions() const
{
  return numOfSubexpress;
}



void TigrRegex::compile()
{
  patbuf.allocated=0;
  patbuf.fastmap=0;
  patbuf.buffer=NULL;
  patbuf.translate=NULL;
  patbuf.syntax=re_syntax_options=
    RE_NO_BK_BRACES  | RE_CHAR_CLASSES |
    RE_NO_BK_PARENS | RE_INTERVALS | 
    RE_NO_BK_VBAR | RE_SYNTAX_EMACS;
  const char *errMsg=re_compile_pattern(regex.c_str(),regex.size(),&patbuf);
  if(errMsg!=NULL) throw TigrString(errMsg);
  numOfSubexpress=patbuf.re_nsub;
}



void TigrRegex::getMatchIndices(int regNum,int &start,int &end)
{
  start=reg.start[regNum];
  end=reg.end[regNum];
}
