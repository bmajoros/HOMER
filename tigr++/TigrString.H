/**************************************************************
TigrString.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrString_H
#define INCL_TigrString_H
#include <stdlib.h>
#include <iostream>
#include <string>
#include "TigrVector.H"
#undef toupper
#undef tolower
using namespace std;

class TigrString : public string
{
public:
  TigrString(unsigned);
  TigrString(unsigned long);
  TigrString(long);
  TigrString(int);
  TigrString(float);
  TigrString(double);
  TigrString(const char *,unsigned len);
  TigrString(const char *);
  TigrString(const TigrString &);
  TigrString(const string &);
  TigrString(char);
  TigrString();
  virtual ~TigrString();
  TigrString operator+(unsigned);
  TigrString operator+(unsigned long);
  TigrString operator+(long);
  TigrString operator+(int);
  TigrString operator+(float);
  TigrString operator+(double);
  TigrString operator+(const char *p);
  TigrString operator+(const string &s);
  TigrString operator+(const TigrString &s);
  TigrString operator+(char c);
  operator float() {return asFloat();}
  operator double() {return asDouble();}
  operator int() {return asInt();}
  operator long() {return asLong();}
  operator unsigned() {return asUnsigned();}
  TigrString &tolower();
  TigrString &toupper();
  TigrString substitute(const TigrString &changeThis,
			const TigrString &toThis) const;
  TigrString substring(int start,int length) const;
  TigrVector<TigrString> *getFields(const char *separators=" \t\n\r") const;
  bool contains(const TigrString &searchString) const;
  bool containsWordChar();
  inline bool isWordChar(char c);
  bool stricmp(const TigrString &) const;
  char lastChar() const { return (*this)[length()-1]; }
  double asDouble() { return atof(c_str()); }
  float asFloat() { return atof(c_str()); }
  int asInt() { return atoi(c_str()); }
  long asLong() { return atol(c_str()); }
  unsigned asUnsigned() { return asLong(); }
  void chop();
  void trimWhitespace();
  inline bool isEmpty() {return empty();}
  void getline(istream &is) {std::getline(is,*this);}
private:
  bool localMatch(const char *,const char *,int) const;
};

ostream &operator<<(ostream &,const TigrString &);


//This is bad!!! Do not uncomment this line!!!
//inline istream &operator>>(istream &is,TigrString &s) {return std::getline(is,s);}


#endif
