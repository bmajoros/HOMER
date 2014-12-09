/**************************************************************
TigrStringMap.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrStringMap_H
#define INCL_TigrStringMap_H
#include <list>
#include "TigrString.H"
using namespace std;

/*****************************************************************
  Some primes to use for table sizes:

  tens:
    13 19 23 31 41 53 61 71 83 89

  hundreds:
    101 199 293 401 499 601 701 797 887 

  thousands:
    997 1097 1201 1301 1399 1499 1601 1699 1801 1901 1999 2099 2179 
    2297 2399 2477 2593 2699 2801 2897 3001 3089 3191 3301 3391 3499 
    3593 3701 3797 3889 4001 4099 4201 4297 4397 4493 4597 4691 4801 
    4889 4999 5101 5197 5297 5399 5501 5591 5701 5801 5897 5987 6101 
    6199 6301 6397 6491 6599 6701 6793 6899 7001 7079 7193 7297 7393 
    7499 7591 7699 7793 7901 7993 8101 8191 8297 8389 8501 8599 8699 
    8783 8893 9001 9091 9199 9293 9397 9497 9601 9697 9791 9901 

  ten-thousands:
    9973 10993 11987 12983 13999 14983 15991 16993 17989 18979 19997
    20983 21997 22993 23993 24989 25999 26993 27997 28979 29989 30983
    31991 32999 33997 34981 35999 36997 37997 38993 39989 40993 41999
    42989 43997 44987 45989 46997 47981 48991 49999 50993 51991 52999
    53993 54983 55997 56999 57991 58997 59999 60961 61991 62989 63997
    64997 65993 66977 67993 68993 69997 70999 71999 72997 73999 74959
    75997 76991 77999 78989 79999 80989 81973 82997 83987 84991 85999
    86993 87991 88997 89989 90997 91997 92993 93997 94999 95989 96997 
    97987 98999 

  hundred-thousands:
    99961 199999 299993 399989 499979 599999 700001 799999 900001 

  millions:
    999959 1999957 
 *****************************************************************/

/*****************************************************************
 */
template<class T>
struct StringMapElem
{
  int len;
  char *first;
  T second;
  StringMapElem(const char *,int len,const T &);
  StringMapElem(const char *,int len);
  StringMapElem(const StringMapElem<T> &);
  virtual ~StringMapElem();
};

/*****************************************************************
 */
template<class T>
class StringMapIterator
{
public:
  typedef list<StringMapElem<T>*> Bucket;
  typedef StringMapElem<T> ElemType;
  StringMapIterator(int index,int tableSize,Bucket *array);
  inline bool operator!=(const StringMapIterator<T> &);
  bool operator==(const StringMapIterator<T> &);
  StringMapIterator<T> &operator++(int);
  StringMapIterator<T> &operator++();
  StringMapElem<T> &operator*();
private:
  int index, tableSize;
  typename Bucket::iterator cur, end;
  Bucket *array;
  void findNonemptyBucket();
};

/*****************************************************************
 */
template<class T>
class TigrStringMap
{
public:
  typedef StringMapElem<T> ElemType;
  typedef StringMapIterator<T> iterator;
  typedef StringMapIterator<T> const_iterator;
  TigrStringMap(int tableSize);
  TigrStringMap(const TigrStringMap<T> &);
  virtual ~TigrStringMap();
  T &lookup(const char *,int index,int len);
  T &lookup(const char *,int len);
  bool isDefined(const char *,int index,int len);
  bool isDefined(const char *,int len);
  int size();
  iterator begin();
  iterator end();
  iterator begin() const;
  iterator end() const;
  void clear();
  void remove(const char *,int index,int len);
  void remove(const char *,int len);
  TigrStringMap<T> &operator=(const TigrStringMap &);
private:
  typedef list<StringMapElem<T>*> Bucket;
  int tableSize, numElements;
  Bucket *array;
  unsigned hash(const char *,int len,int tableSize);
};



/*****************************************************************
 */
template<class T>
StringMapElem<T>::StringMapElem(const StringMapElem<T> &other)
  : len(other.len), second(other.second)
{
  first=new char[len+1];
  memcpy(first,other.first,len+1);
}



template<class T>
StringMapElem<T>::StringMapElem(const char *p,int len)
  : len(len)
{
  first=new char[len+1];
  first[len]='\0';
  strncpy(first,p,len);
}



template<class T>
StringMapElem<T>::StringMapElem(const char *p,int len,const T &t)
  : second(t), len(len)
{
  first=new char[len+1];
  first[len]='\0';
  strncpy(first,p,len);
}



template<class T>
StringMapElem<T>::~StringMapElem()
{
  delete [] first;
}



template<class T>
StringMapElem<T> &StringMapIterator<T>::operator*()
{
  return **cur;  
}



template<class T>
StringMapIterator<T> &StringMapIterator<T>::operator++()
{
  if(index<tableSize)
    {
      ++cur;
      if(cur==end) 
	{
	  ++index;
	  findNonemptyBucket();
	}
    }
  return *this;
}



template<class T>
StringMapIterator<T> &StringMapIterator<T>::operator++(int)
{
  if(index<tableSize)
    {
      ++cur;
      if(cur==end) 
	{
	  ++index;
	  findNonemptyBucket();
	}
    }
  return *this;
}



template<class T>
StringMapIterator<T>::StringMapIterator(int index,int tableSize,Bucket *array)
  : array(array), tableSize(tableSize), index(index)
{
  findNonemptyBucket();  
}



template<class T>
bool StringMapIterator<T>::operator!=(const StringMapIterator<T> &i) 
{
  return !(*this==i);
}



template<class T>
bool StringMapIterator<T>::operator==(const StringMapIterator<T> &other)
{
  if(other.index!=index) return false;
  if(index>=tableSize) return true;
  return cur==other.cur;
}



template<class T>
void StringMapIterator<T>::findNonemptyBucket()
{
  for(; index<tableSize ; ++index)
    {
      Bucket &bucket=array[index];
      if(!bucket.empty())
	{
	  cur=bucket.begin();
	  end=bucket.end();
	  return;
	}
    }
}



template<class T>
T &TigrStringMap<T>::lookup(const char *p,int len)
{
  return lookup(p,0,len);
}



template<class T>
T &TigrStringMap<T>::lookup(const char *pOrigin,int index,int len)
{
  const char *p=pOrigin+index;
  unsigned hashValue=hash(p,len,tableSize);
  Bucket &bucket=array[hashValue];
  typename Bucket::iterator end=bucket.end();
  typename Bucket::iterator cur=bucket.begin();
  for(; cur!=end ; ++cur)
    {
      StringMapElem<T> *elem=*cur;
      if(len!=elem->len) continue;
      if(!strncmp(p,elem->first,len))
	return elem->second;
    }
  StringMapElem<T> *newElem=new StringMapElem<T>(p,len);
  bucket.push_back(newElem);
  ++numElements;
  return newElem->second;
}



template<class T>
TigrStringMap<T>::TigrStringMap(const TigrStringMap<T> &other)
  : tableSize(other.tableSize), numElements(other.numElements),
    array(new Bucket[other.tableSize])
{
  for(int i=0 ; i<tableSize ; ++i)
    {
      Bucket &thisBucket=array[i], &thatBucket=other.array[i];
      typename Bucket::iterator cur=thatBucket.begin(), end=thatBucket.end();
      for(; cur!=end ; ++cur)
	thisBucket.push_back(new StringMapElem<T>(**cur));
    }
}



template<class T>
TigrStringMap<T> &TigrStringMap<T>::operator=(const TigrStringMap &other)
{
  tableSize=other.tableSize;
  numElements=other.numElements;
  array=new Bucket[other.tableSize];
  for(int i=0 ; i<tableSize ; ++i)
    {
      Bucket &thisBucket=array[i];
      const Bucket &thatBucket=other.array[i];
      typename Bucket::const_iterator cur=thatBucket.begin(), 
	end=thatBucket.end();
      for(; cur!=end ; ++cur)
	thisBucket.push_back(new StringMapElem<T>(**cur));
    }
  return *this;
}



template<class T>
TigrStringMap<T>::TigrStringMap(int tableSize)
  : tableSize(tableSize), numElements(0)
{
  array=new Bucket[tableSize];
}



template<class T>
typename TigrStringMap<T>::iterator TigrStringMap<T>::begin()
{
  return iterator(0,tableSize,array);
}



template<class T>
typename TigrStringMap<T>::iterator TigrStringMap<T>::end()
{
  return iterator(tableSize,tableSize,array);
}



template<class T>
typename TigrStringMap<T>::iterator TigrStringMap<T>::begin() const
{
  return iterator(0,tableSize,const_cast<Bucket*>(array));
}



template<class T>
typename TigrStringMap<T>::iterator TigrStringMap<T>::end() const
{
  return iterator(tableSize,tableSize,const_cast<Bucket*>(array));
}



template<class T>
TigrStringMap<T>::~TigrStringMap()
{
  //cout<<"~TigrStringMap"<<endl;
  for(int i=0 ; i<tableSize ; ++i)
    {
      Bucket &bucket=array[i];
      typename Bucket::iterator cur=bucket.begin(), end=bucket.end();
      for(; cur!=end ; ++cur)
	delete *cur;
    }
  delete [] array;
}



template<class T>
bool TigrStringMap<T>::isDefined(const char *p,int len)
{
  return isDefined(p,0,len);
}



template<class T>
bool TigrStringMap<T>::isDefined(const char *pOrigin,int index,int len)
{
  const char *p=pOrigin+index;
  unsigned hashValue=hash(p,len,tableSize);
  Bucket &bucket=array[hashValue];
  typename Bucket::iterator cur=bucket.begin(), end=bucket.end();
  for(; cur!=end ; ++cur)
    {
      StringMapElem<T> *elem=*cur;
      if(len!=elem->len) continue;
      if(!strncmp(p,elem->first,len))
	return true;
    }
  return false;
}



template<class T>
int TigrStringMap<T>::size()
{
  return numElements;
}



template<class T>
unsigned TigrStringMap<T>::hash(const char *s,int length,int tableSize)
{
  int h=0;
  const char *p=s, *end=s+length;
  for(; p!=end ; ++p)
    {
      h=(h<<4)+*p;
      int g=h & 0xf000;
      if(g) h=h^(g>>8);
    }
  return labs(h) % tableSize;
}



template<class T>
void TigrStringMap<T>::clear()
{
  for(int i=0 ; i<tableSize ; ++i)
    {
      Bucket &bucket=array[i];
      typename Bucket::iterator cur=bucket.begin(), end=bucket.end();
      for(; cur!=end ; ++cur)
	delete *cur;
      bucket.clear();
    }
  numElements=0;
}



template<class T>
void TigrStringMap<T>::remove(const char *p,int len)
{
  remove(p,0,len);
}



template<class T>
void TigrStringMap<T>::remove(const char *pOrigin,int index,int len)
{
  const char *p=pOrigin+index;
  unsigned hashValue=hash(p,len,tableSize);
  Bucket &bucket=array[hashValue];
  typename Bucket::iterator cur=bucket.begin(), end=bucket.end();
  for(; cur!=end ; ++cur)
    {
      StringMapElem<T> *elem=*cur;
      if(len!=elem->len) continue;
      if(!strncmp(p,elem->first,len))
	{
	  bucket.erase(cur);
	  break;
	}
    }
  --numElements;
}

#endif
