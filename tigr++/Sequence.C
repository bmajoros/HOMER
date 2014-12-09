/**************************************************************
Sequence.C
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
using namespace std;
#include "Sequence.H"
#include "TigrFastaReader.H"


TigrRegex Sequence::deflineRegex("^\\s*>\\s*(\\S+).*");


Sequence::Sequence(const TigrString &s,Alphabet &alphabet)
  : phase(-1)
{
  append(alphabet,s.c_str());
}



Sequence::Sequence() : phase(-1)
{
}



Sequence::Sequence(const char *str,Alphabet &alphabet)
  : phase(-1)
{
  append(alphabet,str);
}



void Sequence::copyFrom(const TigrString &source,Alphabet &alphabet)
{
  phase=-1;
  symbols.clear();
  append(alphabet,source.c_str());
}



void Sequence::append(Alphabet &alphabet,const char *ptr)
{
  while(*ptr)
    if(!isspace(*ptr))
      {
	const char c=*ptr++;
	Symbol s=alphabet.lookup(c);
	append(s);
      }
    else
      ptr++;
}



ostream &operator<<(ostream &os,const Sequence &s)
{
  Sequence &ncs=const_cast<Sequence&>(s);
  ncs.printOn(os);
  return os;
}



void Sequence::printOn(ostream &os) const
{
  int numSymbols=symbols.size();
  for(int i=0 ; i<numSymbols ; ++i) os<<symbols[i]<<' ';
}



int Sequence::getLength() const
{
  return symbols.size();
}



void Sequence::printOn(ostream &os,Alphabet &alphabet) const
{
  TigrString *s=toString(alphabet);
  os << *s;
  delete s;
}



TigrString *Sequence::toString(Alphabet &alphabet,int startingAt) const
{
  TigrString &str=*new TigrString;
  int numSymbols=symbols.size();
  for(int i=startingAt ; i<numSymbols ; ++i) 
    str+=alphabet.lookup(symbols[i]);
  return &str;
}



void Sequence::clear()
{
  symbols.clear();
}



void Sequence::getSubsequence(int begin,int len,Sequence &seq) const
{
  const Sequence &self=*this;
  seq.clear();
  int last=begin+len;
  for(int i=begin ; i<last ; ++i)
    seq.append(self[i]);
}



bool Sequence::subsequenceOccursAt(const Sequence &subsequence,int at) const
{
  const Sequence &self=*this;
  int end=at+subsequence.getLength();
  if(end>=getLength()) return false;
  for(int i=at, j=0 ; i<end ; ++i, ++j)
    if(self[i]!=subsequence[j])
      return false;
  return true;
}



void Sequence::setPhase(int phase)
{
  this->phase=phase;
}



Symbol &Sequence::operator[](int i)
{
  return symbols[i];
}



void Sequence::prepend(Symbol s)
{
  symbols.push_front(s);
}



void Sequence::append(Symbol s)
{
  symbols.push_back(s);
}



void Sequence::append(const Sequence &other)
{
  symbols.append(other.symbols);
}



int Sequence::getPhase() const
{
  return phase;
}



Symbol Sequence::operator [](int i) const
{
  return symbols[i];
}



Sequence *Sequence::load(const TigrString &filename,Alphabet &alphabet,
			 TigrString &substrateId)
{
  TigrFastaReader reader(filename,alphabet);
  TigrString sequence, defline;
  reader.nextSequence(defline,sequence);
  if(!deflineRegex.match(defline))
    throw TigrString("Can't parse defline: \"")+defline+"\"";
  substrateId=deflineRegex[1];
  return new Sequence(sequence,alphabet);
}



Sequence *Sequence::load(const TigrString &filename,Alphabet &alphabet)
{
  TigrString substrateId;
  return load(filename,alphabet,substrateId);
}



const Sequence &Sequence::operator=(const Sequence &s)
{
  clear();
  append(s);
  return s;
}




void Sequence::reverseComplement(Alphabet &alphabet,Sequence &rc) const
{
  rc.phase=(phase+getLength()-1) % 3;

  TigrVector<Symbol>::const_reverse_iterator cur=symbols.rbegin(), 
    end=symbols.rend();
  for(; cur!=end ; ++cur)
    rc.append(alphabet.complement(*cur));
}



Sequence *Sequence::reverseComplement(Alphabet &alphabet) const
{
  Sequence *rc=new Sequence;
  reverseComplement(alphabet,*rc);
  return rc;
}



int Sequence::countOccurrences(Symbol s)
{
  int c=0, n=symbols.size();
  for(int i=0 ; i<n ; ++i)
    if(symbols[i]==s)
      ++c;
  return c;
}

