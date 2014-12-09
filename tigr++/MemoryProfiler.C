/****************************************************************
 MemoryProfiler.C
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
using namespace std;
#include "MemoryProfiler.H"
#include <iostream>
#include "unistd.h"
#include "TigrFile.H"


void MemoryProfiler::report(const TigrString &label,ostream &os)
{
  int pid=getpid();
  TigrString filename=TigrString("/proc/")+pid+TigrString("/status");
  TigrFile file(filename,"r");
  long size=-1;
  TigrString units;
  while(!file.eof())
    {
      TigrString line=file.getline();
      if(file.eof()) break;
      TigrVector<TigrString> &fields=*line.getFields();
      if(fields.size()>0 && fields[0]=="VmSize:")
	{
	  size=fields[1].asLong();
	  units=fields[2];
	  delete &fields;
	  break;
	}
      delete &fields;
    }
  os<<label<<" ";
  if(units=="kB" && size>=1000)
    os<<(size/1000)<<" Mb"<<endl;  
  else 
    os<<size<<" "<<units<<endl;
}


