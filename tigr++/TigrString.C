#include <sstream>
//#include <strstream>
#include <iostream>
#include "TigrString.H"
#include "TigrStrTokenizer.H"
using namespace std;



TigrString::TigrString(unsigned u) : string(TigrString("")+u)
{
}



TigrString::TigrString(unsigned long ul) : string(TigrString("")+ul)
{
}



TigrString::TigrString(long l) : string(TigrString("")+l)
{
}



TigrString::TigrString(int i) : string(TigrString("")+i)
{
}



TigrString::TigrString(float f) : string(TigrString("")+f)
{
}



TigrString::TigrString(double d) : string(TigrString("")+d)
{
}



TigrString::TigrString(const char *cp,unsigned length) : string(cp,length)
{
}



TigrString::TigrString(const char *p) : string(p?p:"")
{
}



TigrString::TigrString(const TigrString &s) : string(s)
{
}



TigrString::TigrString(const string &s) : string(s)
{
}



TigrString::TigrString(char c) : string(TigrString("")+c)
{
}



TigrString::TigrString() {}


TigrString::~TigrString() {}


ostream &operator<<(ostream &os,const TigrString &str)
{
  os << str.c_str();
  return os;
}



TigrString TigrString::operator+(unsigned u)
{
  stringstream os;
  os << c_str() << u;
  return os.str();
}



TigrString TigrString::operator+(unsigned long ul)
{
  stringstream os;
  os << c_str() << ul;
  return os.str();
}



TigrString TigrString::operator+(long l)
{
  stringstream os;
  os << c_str() << l;
  return os.str();
}



TigrString TigrString::operator+(int i)
{
  stringstream os;
  os << c_str() << i;
  return os.str();
}



TigrString TigrString::operator+(float f)
{
  stringstream os;
  os << c_str() << f;
  return os.str();
}



TigrString TigrString::operator+(double d)
{
  stringstream os;
  os << c_str() << d;
  return os.str();
}



TigrString &TigrString::tolower()
{
  TigrString &self=*this;
  int l=length();
  for(int i=0 ; i<l ; ++i)
    self[i]=::tolower(self[i]);
  return *this;
}



TigrString &TigrString::toupper()
{
  TigrString &self=*this;
  int len=length();
  for(int i=0 ; i<len ; ++i)
    self[i]=::toupper(self[i]);
  return *this;
}



TigrString TigrString::substitute(const TigrString &from,
				  const TigrString &to) const
{
  TigrString rval;
  const char *pattern=from.c_str();
  int patternLen=from.length();
  const char *ptr=c_str();
  const char *last=ptr+length()-patternLen;
  while(ptr<=last)
    {
      if(localMatch(ptr,pattern,patternLen))
	{
	  ptr+=patternLen;
	  rval+=to;
	}
      else
	{
	  rval+=*ptr;
	  ptr++;
	}
    }
  return rval;
}



TigrString TigrString::substring(int begin,int len) const
{
  return substr(begin,len);
}



TigrVector<TigrString> *TigrString::getFields(const char *seps) const
{
  TigrString &self=const_cast<TigrString&>(*this);
  TigrStrTokenizer tokenizer(self,seps);
  TigrVector<TigrString> *fields=tokenizer.getTokenStrings();
  return fields;
}



bool TigrString::contains(const TigrString &s) const
{
  return find(s.c_str())!=npos;
}



bool TigrString::isWordChar(char c)
{
  return (c>='A' && c<='Z') || (c>='a' && c<='z');
}



bool TigrString::containsWordChar()
{
  for(const char *ptr=c_str() ; *ptr ; ++ptr)
    {
      char c=*ptr;
      if(isWordChar(c)) return true;
    }
  return false;
}



bool TigrString::localMatch(const char *s1,const char *s2,int len) const
{
  for(int i=0 ; i<len ; ++i)
    if(*s1++!=*s2++) 
      return false;
  return true;
}



bool TigrString::stricmp(const TigrString &str) const
{
  return strcasecmp(c_str(),str.c_str());
}



void TigrString::chop()
{
  int len=length();
  resize(len-1);
}



void TigrString::trimWhitespace()
{
  TigrString &self=*this;
  const char *begin=c_str();
  while(*begin && isspace(*begin)) 
    ++begin;
  const char *end=begin;
  while(*end) 
    ++end;
  if(end>begin) 
    --end;
  while(end>begin && isspace(*end)) 
    --end;
  int len=end-begin+2;
  char *buf=new char[len];
  char *ptr=buf;
  while(begin<=end) 
    *ptr++=*begin++;
  *ptr='\0';
  self=buf;
  delete [] buf;
}



TigrString TigrString::operator+(const char *p) 
{
  return TigrString(((string&)*this)+p);
}



TigrString TigrString::operator+(const string &s) 
{
  return TigrString(*this+s.c_str());
}



TigrString TigrString::operator+(const TigrString &s) 
{
  return TigrString(*this+s.c_str());
}



TigrString TigrString::operator+(char c)
{
  char str[2];
  str[0]=c;
  str[1]='\0';
  return TigrString(*this+(const char*)str);
}
