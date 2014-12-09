/**************************************************************
TigrFile.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrFile_H
#define INCL_TigrFile_H
#include <unistd.h>
#include <sys/stat.h>
#include <stdio.h>
#include "TigrString.H"
using namespace std;

class TigrFile
{
public:
  TigrFile(const TigrString &filename,const TigrString &mode="r");
  TigrFile(FILE *fp);
  TigrFile();
  virtual ~TigrFile();
  operator FILE*();
  long countLines();
  static TigrString getFilenameNoPath(const TigrString &filename);
  static bool copy(const TigrString &from,const TigrString &to);
  static bool exists(const TigrString &filename);
  TigrString getline();  
  inline TigrString readLine() { return getline(); }
  bool eof();
  bool isOpen() const;
  bool open();
  bool open(const TigrString &filename,const TigrString &mode="r");
  long read(long bytes,void *buffer);
  bool seek(long offset);
  bool write(long bytes,const void *buffer);
  long getPosition();
  long getSize();
  time_t lastAccessTime();
  time_t lastChangeTime();
  time_t lastModifyTime();
  virtual void close();
  void print(const TigrString &);
  void rewind();
  void updateStats();
  void write(const int &);
  void write(const short &);
  void write(const long &);
  void write(char);
  void write(const float &);
  void write(const double &);
  void write(const TigrString &);
  void write(const char *);
  int readInt();
  short readShort();
  long readLong();
  char readChar();
  float readFloat();
  double readDouble();
  TigrString readString();
  char *readCharString();
  void readString(TigrString &);
  TigrString *readStringPtr();
protected:
  struct stat fileInfo;
  TigrString filename, mode;
  FILE *fp;
  static const char *findLastSlash(const char *begin,const char *end);
  static const char *findEndOfString(const char *begin);
};

// Writing:
TigrFile &operator<<(TigrFile &,int);
TigrFile &operator<<(TigrFile &,short);
TigrFile &operator<<(TigrFile &,long);
TigrFile &operator<<(TigrFile &,float);
TigrFile &operator<<(TigrFile &,double);
TigrFile &operator<<(TigrFile &,char);
TigrFile &operator<<(TigrFile &,const TigrString &);
TigrFile &operator<<(TigrFile &,const char *);

// Reading:
TigrFile &operator>>(TigrFile &,int &);
TigrFile &operator>>(TigrFile &,short &);
TigrFile &operator>>(TigrFile &,long &);
TigrFile &operator>>(TigrFile &,float &);
TigrFile &operator>>(TigrFile &,double &);
TigrFile &operator>>(TigrFile &,char &);
TigrFile &operator>>(TigrFile &,TigrString &);

#endif
