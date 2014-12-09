/**************************************************************
TigrBitSet.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrBitSet_H
#define INCL_TigrBitSet_H
using namespace std;
#include <stdio.h>

class TigrBitSet 
{
public:
  TigrBitSet(unsigned long Size=0);  
  virtual ~TigrBitSet();
  void operator =(TigrBitSet &RValue);   
  void operator -=(TigrBitSet &otherSet);  
  void operator +=(TigrBitSet &otherSet);  
  void operator *=(TigrBitSet &otherSet);  
  bool operator ==(TigrBitSet &otherSet);   
  TigrBitSet *operator -(TigrBitSet &otherSet);  
  TigrBitSet *operator +(TigrBitSet &otherSet);  
  TigrBitSet *operator *(TigrBitSet &otherSet);  
  void intersect(TigrBitSet &otherSet,TigrBitSet &into);
  bool isMember(unsigned long BitNumber);  
  unsigned long cardinality();
  unsigned long getMaxSize();
  void addAll();  
  void addMember(unsigned long BitNumber);  
  void complement();
  void getRawBytes(unsigned char *&,unsigned long &maxSize,
		   unsigned long &numBytes);
  void load(FILE *fp);
  void purge();                
  void removeMember(unsigned long BitNumber);  
  void replaceRawBytes(unsigned char *,unsigned long maxSize,
		       unsigned long numBytes);
  void save(FILE *fp);  
  void setSize(unsigned long Size);
protected:
  unsigned long NumBytes;  
  unsigned long MaxSize;   
  unsigned char *TheBitMap;
};

#endif
