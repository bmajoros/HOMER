using namespace std;
#include "TigrExceptions.H"
#include "TigrBitSet.H"



TigrBitSet::TigrBitSet(unsigned long Size) 
  : TheBitMap(NULL)
{
  setSize(Size);
}



TigrBitSet::~TigrBitSet() 
{ 
  delete [] TheBitMap; 
}



void TigrBitSet::operator =(TigrBitSet &RValue) 
{  
  MaxSize=RValue.MaxSize;
  NumBytes=RValue.NumBytes;
  
  delete [] TheBitMap;
  TheBitMap=new unsigned char[NumBytes];

  unsigned long f;
  for(f=0 ; f<NumBytes ; f++)
    TheBitMap[f]=RValue.TheBitMap[f];
}



void TigrBitSet::operator -=(TigrBitSet &OtherSet) 
{
  unsigned long f;
  for(f=NumBytes ; f>0 ; f--)
    TheBitMap[f-1] &= (~(OtherSet.TheBitMap[f]-1));
}



void TigrBitSet::operator +=(TigrBitSet &OtherSet) 
{
  unsigned long f;
  for(f=NumBytes ; f>0 ; f--)
    TheBitMap[f-1] |= OtherSet.TheBitMap[f-1];
}



void TigrBitSet::operator *=(TigrBitSet &OtherSet) 
{
  unsigned long f;
  for(f=NumBytes ; f>0 ; f--)
    TheBitMap[f-1] &= OtherSet.TheBitMap[f-1];
}



bool TigrBitSet::operator ==(TigrBitSet &OtherSet) 
{
  unsigned long f;
  for(f=NumBytes-1 ; f ; f--)
    if(TheBitMap[f]!=OtherSet.TheBitMap[f]) 
      return false;
  
  unsigned char BitMask=1;
  unsigned long RelevantBits=MaxSize % 8;
  for(f=1; f<=RelevantBits; f++) 
    {
      if((BitMask & TheBitMap[0])!=
	 (BitMask & OtherSet.TheBitMap[0])) 
	return false;
      BitMask <<=1;
    }
  
  return true;
}



TigrBitSet *TigrBitSet::operator -(TigrBitSet &OtherSet) 
{
  TigrBitSet *retval=new TigrBitSet(MaxSize);
  *retval=*this;
  *retval-=OtherSet;
  
  return retval;
}



TigrBitSet *TigrBitSet::operator +(TigrBitSet &OtherSet) 
{
  TigrBitSet *retval=new TigrBitSet(MaxSize);
  *retval+=OtherSet;
  *retval+=*this;
  
  return retval;
}



void TigrBitSet::intersect(TigrBitSet &otherSet,TigrBitSet &retval)
{
  unsigned char *r=retval.TheBitMap, *t=TheBitMap, *o=otherSet.TheBitMap;
  for(unsigned long f=0 ; f<NumBytes ; ++f)
    *r++ = *t++ & *o++;
}



TigrBitSet *TigrBitSet::operator *(TigrBitSet &otherSet) 
{
  TigrBitSet *retval=new TigrBitSet(MaxSize);
  intersect(otherSet,*retval);
  return retval;
}



bool TigrBitSet::isMember(unsigned long BitNumber) 
{
  if(BitNumber>=MaxSize) 
    throw ArrayIndexException(BitNumber,
       "in TigrBitSet::IsMember()");
  
  unsigned long ByteFromRight=BitNumber/8;
  unsigned long BitInByte=BitNumber%8;

  unsigned char ByteMask=1;
  ByteMask <<= BitInByte;
  
  return (TheBitMap[NumBytes-1-ByteFromRight] & ByteMask) ? 
    true : false;
}



unsigned long TigrBitSet::cardinality() 
{
  unsigned long TotalElements=0, f;
  for(f=NumBytes-1 ; f ; f--) 
    {
      unsigned char BitMask=1;
      unsigned long ElementsInByte=0;
      unsigned long g;
      for(g=1 ; g<=8 ; g++) 
	{
	  if(BitMask & TheBitMap[f]) ElementsInByte++;
	  BitMask <<=1;
	}
      TotalElements+=ElementsInByte;
    }
   
  unsigned long RelevantBits=MaxSize % 8;
  unsigned char BitMask=1;
  unsigned long ElementsInByte=0;
  for(f=1 ; f<=RelevantBits ; f++) 
    {
      if(BitMask & TheBitMap[0]) ElementsInByte++;
      BitMask <<=1;
    }
  TotalElements+=ElementsInByte;
  
  return TotalElements;
}



unsigned long TigrBitSet::getMaxSize() 
{ 
  return MaxSize; 
}



void TigrBitSet::addAll()
{
  unsigned long i;
  for(i=0 ; i<NumBytes ; i++)
    TheBitMap[i]=char(0xFF);  
}



void TigrBitSet::addMember(unsigned long BitNumber) 
{
  if(BitNumber>=MaxSize) 
    throw ArrayIndexException(BitNumber,
			      "in TigrBitSet::AddMember()");
  
  unsigned long ByteFromRight=BitNumber/8;
  unsigned long BitInByte=BitNumber%8;

  unsigned char ByteMask=1;
  ByteMask <<= BitInByte;
  
  TheBitMap[NumBytes-1-ByteFromRight] =
    TheBitMap[NumBytes-1-ByteFromRight] | ByteMask;
}



void TigrBitSet::complement()
{
  unsigned long i;
  for(i=0 ; i<NumBytes ; i++)
    TheBitMap[i]=~TheBitMap[i];  
}



void TigrBitSet::getRawBytes(unsigned char *&bytes,
			    unsigned long &maxSize,
			    unsigned long &numBytes)
{
  bytes=TheBitMap;
  maxSize=this->MaxSize;
  numBytes=this->NumBytes;
}



void TigrBitSet::load(FILE *fp)
{
   
  if(fread(&MaxSize,1,sizeof(MaxSize),fp) != sizeof(MaxSize)) 
    throw FileErrorException("<filename unknown>",
			     "Cannot load TigrBitSet");
  setSize(MaxSize);
  
   
  if(fread(TheBitMap,1,NumBytes,fp) != NumBytes) 
    throw FileErrorException("<filename unknown>",
			     "Cannot load TigrBitSet");
}



void TigrBitSet::purge() 
{
  unsigned long i;
  
  for(i=0 ; i<NumBytes ; i++)
    TheBitMap[i]='\0';
}



void TigrBitSet::removeMember(unsigned long BitNumber) 
{
  if(BitNumber>=MaxSize) 
    throw ArrayIndexException(BitNumber,
			      "in TigrBitSet::RemoveMember()");
  
  unsigned long ByteFromRight=BitNumber/8;
  unsigned long BitInByte=BitNumber%8;

  unsigned char ByteMask=1;
  ByteMask <<= BitInByte;
  
  TheBitMap[NumBytes-1-ByteFromRight] =
    TheBitMap[NumBytes-1-ByteFromRight] & (~ByteMask);
}



void TigrBitSet::replaceRawBytes(unsigned char *rawBytes,
				unsigned long maxSize,
				unsigned long numBytes)
{
  delete [] this->TheBitMap;

  this->TheBitMap=rawBytes;
  this->MaxSize=maxSize;
  this->NumBytes=numBytes;
}



void TigrBitSet::save(FILE *fp)
{
   
  if(fwrite(&MaxSize,1,sizeof(MaxSize),fp)!=sizeof(MaxSize)) 
    throw FileErrorException("<filename unknown>",
			     "Cannot save TigrBitSet");
  
  if(fwrite(TheBitMap,1,NumBytes,fp) != NumBytes) 
    throw FileErrorException("<filename unknown>",
			     "Cannot save TigrBitSet");
}



void TigrBitSet::setSize(unsigned long Size) 
{
  delete [] TheBitMap;
  
  MaxSize=Size;
  NumBytes=MaxSize/8;
  if(MaxSize%8) NumBytes++;
  
  TheBitMap=new unsigned char[NumBytes];
  
  purge();
}
