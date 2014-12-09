#include <iostream>
#include "TigrFile.H"
using namespace std;



TigrFile::TigrFile(const TigrString &filename,const TigrString &mode)
  : mode(mode), filename(filename), fp(NULL) 
{
  if(mode.length()>0) 
    if(mode[0]=='r') 
      updateStats();

  if(!open()) throw TigrString("Cannot open file: \"")+filename+"\"";
}



TigrFile::TigrFile(FILE *fp)
  : fp(fp)
{
}



TigrFile::TigrFile()
  : fp(NULL)
{
}



TigrFile::~TigrFile()
{
  close();
}



TigrString TigrFile::getFilenameNoPath(const TigrString &filename)
{
  const char *begin=filename.c_str();
  const char *end=findEndOfString(begin);
  const char *slash=findLastSlash(begin,end);
  return TigrString(slash);
}



const char *TigrFile::findLastSlash(const char *begin,const char *end)
{
  const char *slash=end;
  while(slash>begin && *slash!='/') 
    --slash;
  if(*slash=='/') 
    ++slash;
  return slash;
}



const char *TigrFile::findEndOfString(const char *end)
{
  while(*end) ++end;
  return end;
}



TigrString TigrFile::getline()
{
  char buf[1025];
  TigrString s;
  while(1)
    {
      if(!fgets(buf,sizeof(buf),fp)) break;
      size_t len=strlen(buf);
      if(buf[len-1]=='\n')
	{
	  buf[len-1]='\0';
	  s+=buf;
	  break;
	}
      s+=buf;
    }
  return s;
}



bool TigrFile::copy(const TigrString &from,const TigrString &to)
{
  TigrString command="cp "+from+" "+to;
  int rval=system(command.c_str())==0;
  return rval;
}



bool TigrFile::eof()
{
  return feof(fp);
}



bool TigrFile::exists(const TigrString &filename)
{
  struct stat fInfo;
  int status=stat(filename.c_str(),&fInfo);
  bool doesExist=(status==0);
  return doesExist;
}



bool TigrFile::isOpen() const
{
  return (bool)(int)fp;
}



bool TigrFile::open()
{
  if(mode.length()>0)
    if(mode[0]=='r') 
      stat(filename.c_str(),&fileInfo);

  if(!fp)
    fp=fopen(filename.c_str(),mode.c_str());
  return (bool)(int)fp;
}



bool TigrFile::open(const TigrString &fn,const TigrString &md)
{
  close();
  if(md.length()>0)
    if(md[0]=='r') 
      stat(fn.c_str(),&fileInfo);
  filename=fn;
  mode=md;
  return open();
}



long TigrFile::read(long numbytes,void *buffer)
{
  return fread(buffer,1,numbytes,fp);
}



bool TigrFile::seek(long pos)
{
  return fseek(fp,pos,SEEK_SET)==0;
}



bool TigrFile::write(long numbytes,const void *buffer)
{
  return fwrite(buffer,1,numbytes,fp)==numbytes;
}



long TigrFile::countLines()
{
  long numlines=0;
  while(!eof())
    {
      getline();
      ++numlines;
    }
  return numlines;
}



long TigrFile::getPosition()
{
  return ftell(fp);
}



long TigrFile::getSize()
{
  return (long) fileInfo.st_size;
}



time_t TigrFile::lastAccessTime()
{
  return fileInfo.st_atime;
}



time_t TigrFile::lastChangeTime()
{
  return fileInfo.st_ctime;
}



time_t TigrFile::lastModifyTime()
{
  return fileInfo.st_mtime;
}



void TigrFile::close()
{
  if(fp) fclose(fp);
  fp=(FILE*)0;
}



void TigrFile::print(const TigrString &s)
{
  if(!fputs(s.c_str(),fp))
    throw TigrString("Error printing into file in TigrFile::print()");
}



void TigrFile::rewind()
{
  ::rewind(fp);
}



void TigrFile::updateStats()
{
  stat(filename.c_str(),&fileInfo);
}



inline TigrFile::operator FILE*() 
{
  return fp;
}



void TigrFile::write(const int &x)
{
  write(sizeof(x),static_cast<const void*>(&x));
}



void TigrFile::write(const short &x)
{
  write(sizeof(x),static_cast<const void*>(&x));
}



void TigrFile::write(const long &x)
{
  write(sizeof(x),static_cast<const void*>(&x));
}



void TigrFile::write(char c)
{
  write(sizeof(c),static_cast<const void*>(&c));
}



void TigrFile::write(const float &x)
{
  write(sizeof(x),static_cast<const void*>(&x));
}



void TigrFile::write(const double &x)
{
  write(sizeof(x),static_cast<const void*>(&x));
}



void TigrFile::write(const TigrString &x)
{
  long len=x.length();
  write(len);
  write(len,x.c_str());
}



void TigrFile::write(const char *p)
{
  long len=strlen(p);
  write(len);
  write(len,p);
}



int TigrFile::readInt()
{
  int x;
  read(sizeof(x),static_cast<void*>(&x));
  return x;
}



short TigrFile::readShort()
{
  short x;
  read(sizeof(x),static_cast<void*>(&x));
  return x;
}



long TigrFile::readLong()
{
  long x;
  read(sizeof(x),static_cast<void*>(&x));
  return x;
}



char TigrFile::readChar()
{
  char x;
  read(sizeof(x),static_cast<void*>(&x));
  return x;
}



float TigrFile::readFloat()
{
  float x;
  read(sizeof(x),static_cast<void*>(&x));
  return x;
}



double TigrFile::readDouble()
{
  double x;
  read(sizeof(x),static_cast<void*>(&x));
  return x;
}



TigrString TigrFile::readString()
{
  int len=readLong();
  char *p=new char[len+1];
  p[len]='\0';
  read(len,p);
  TigrString s(p);
  delete [] p;
  return s;
}



char *TigrFile::readCharString()
{
  int len=readLong();
  char *p=new char[len+1];
  p[len]='\0';
  read(len,p);
  return p;
}



void TigrFile::readString(TigrString &str)
{
  int len=readLong();
  char *p=new char[len+1];
  p[len]='\0';
  read(len,p);
  str=p;
  delete [] p;
}



TigrString *TigrFile::readStringPtr()
{
  int len=readLong();
  char *p=new char[len+1];
  p[len]='\0';
  read(len,p);
  TigrString *str=new TigrString(p);
  delete [] p;
  return str;
}



TigrFile &operator<<(TigrFile &f,int x)
{
  f.write(x);
  return f;
}



TigrFile &operator<<(TigrFile &f,short x)
{
  f.write(x);
  return f;
}



TigrFile &operator<<(TigrFile &f,long x)
{
  f.write(x);
  return f;
}



TigrFile &operator<<(TigrFile &f,float x)
{
  f.write(x);
  return f;
}



TigrFile &operator<<(TigrFile &f,double x)
{
  f.write(x);
  return f;
}



TigrFile &operator<<(TigrFile &f,char x)
{
  f.write(x);
  return f;
}



TigrFile &operator<<(TigrFile &f,const TigrString &x)
{
  f.write(x);
  return f;
}



TigrFile &operator<<(TigrFile &f,const char *x)
{
  f.write(x);
  return f;
}



TigrFile &operator>>(TigrFile &f,int &x)
{
  x=f.readInt();
  return f;
}



TigrFile &operator>>(TigrFile &f,short &x)
{
  x=f.readShort();
  return f;
}



TigrFile &operator>>(TigrFile &f,long &x)
{
  x=f.readLong();
  return f;
}



TigrFile &operator>>(TigrFile &f,float &x)
{
  x=f.readFloat();
  return f;
}



TigrFile &operator>>(TigrFile &f,double &x)
{
  x=f.readDouble();
  return f;
}



TigrFile &operator>>(TigrFile &f,char &x)
{
  x=f.readChar();
  return f;
}



TigrFile &operator>>(TigrFile &f,TigrString &x)
{
  f.readString(x);
  return f;
}


