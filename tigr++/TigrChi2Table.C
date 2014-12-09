#include <iostream>
#include <fstream>
#include <sstream>
#include "TigrVector.H"
#include "TigrFile.H"
#include "TigrChi2Table.H"
using namespace std;

const char *embeddedTable="df 0.25 0.20 0.15 0.10 0.05 0.025 0.02 0.01 0.005 0.0025 0.001 0.0005\n\
1 1.32 1.64 2.07 2.71 3.84 5.02 5.41 6.63 7.88 9.14 10.83 12.12\n\
2 2.77 3.22 3.79 4.61 5.99 7.38 7.82 9.21 10.60 11.98 13.82 15.20\n\
3 4.11 4.64 5.32 6.25 7.81 9.35 9.84 11.34 12.84 14.32 16.27 17.73\n\
4 5.39 5.59 6.74 7.78 9.49 11.14 11.67 13.23 14.86 16.42 18.47 20.00\n\
5 6.63 7.29 8.12 9.24 11.07 12.83 13.33 15.09 16.75 18.39 20.51 22.11\n\
6 7.84 8.56 9.45 10.64 12.53 14.45 15.03 16.81 18.55 20.25 22.46 24.10\n\
7 9.04 9.80 10.75 12.02 14.07 16.01 16.62 18.48 20.28 22.04 24.32 26.02\n\
8 10.22 11.03 12.03 13.36 15.51 17.53 18.17 20.09 21.95 23.77 26.12 27.87\n\
9 11.39 12.24 13.29 14.68 16.92 19.02 19.63 21.67 23.59 25.46 27.83 29.67\n\
10 12.55 13.44 14.53 15.99 18.31 20.48 21.16 23.21 25.19 27.11 29.59 31.42\n\
11 13.70 14.63 15.77 17.29 19.68 21.92 22.62 24.72 26.76 28.73 31.26 33.14\n\
12 14.85 15.81 16.99 18.55 21.03 23.34 24.05 26.22 28.30 30.32 32.91 34.82\n\
13 15.98 16.58 18.90 19.81 22.36 24.74 25.47 27.69 29.82 31.88 34.53 36.48\n\
14 17.12 18.15 19.4 21.06 23.68 26.12 26.87 29.14 31.32 33.43 36.12 38.11\n\
15 18.25 19.31 20.60 22.31 25.00 27.49 28.26 30.58 32.80 34.95 37.70 39.72\n\
16 19.37 20.47 21.79 23.54 26.30 28.85 29.63 32.00 34.27 36.46 39.25 41.31\n\
17 20.49 21.61 22.98 24.77 27.59 30.19 31.00 33.41 35.72 37.95 40.79 42.88\n\
18 21.60 22.76 24.16 25.99 28.87 31.53 32.35 34.81 37.16 39.42 42.31 44.43\n\
19 22.72 23.90 25.33 27.20 30.14 32.85 33.69 36.19 38.58 40.88 43.82 45.97\n\
20 23.83 25.04 26.50 28.41 31.41 34.17 35.02 37.57 40.00 42.34 45.31 47.50\n\
21 24.93 26.17 27.66 29.62 32.67 35.48 36.34 38.93 41.40 43.78 46.80 49.01\n\
22 26.04 27.30 28.82 30.81 33.92 36.78 37.66 40.29 42.80 45.20 48.27 50.51\n\
23 27.14 28.43 29.98 32.01 35.17 38.08 38.97 41.64 44.18 46.62 49.73 52.00\n\
24 28.24 29.55 31.13 33.20 36.42 39.36 40.27 42.98 45.56 48.03 51.18 53.48\n\
25 29.34 30.68 32.28 34.38 37.65 40.65 41.57 44.31 46.93 49.44 52.62 54.95\n\
26 30.43 31.79 33.43 35.56 38.89 41.92 42.86 45.64 48.29 50.83 54.05 56.41\n\
27 31.53 32.91 34.57 36.74 40.11 43.19 44.14 46.96 49.64 52.22 55.48 57.86\n\
28 32.62 34.03 35.71 37.92 41.34 44.46 45.42 48.28 50.99 53.59 56.89 59.30\n\
29 33.71 35.14 36.85 39.09 42.56 45.72 46.69 49.59 52.34 54.97 58.30 60.73\n\
30 34.80 36.25 37.99 40.26 43.77 46.98 47.96 50.89 53.67 56.33 59.70 62.16\n\
40 45.62 47.27 49.24 51.81 55.76 59.34 60.44 63.69 66.77 69.70 73.40 76.09\n\
50 56.33 58.16 60.35 63.17 67.50 71.42 72.61 76.15 79.49 82.66 86.66 89.56\n\
60 66.98 68.97 71.34 74.40 79.08 83.30 84.58 88.38 91.95 95.34 99.61 102.7\n\
80 88.13 90.41 93.11 96.58 101.9 106.6 108.1 112.3 116.3 120.1 124.8 128.3\n\
100 109.1 111.7 114.7 118.5 124.3 129.6 131.1 135.8 140.2 144.3 149.4 153.2\n\
";



TigrChi2Table::TigrChi2Table(TigrString filename)
  : DFs(35), alphas(12), table(12,35)
{
  if(filename=="")
    {
      istringstream is(embeddedTable);
      loadFromStream(is);
    }
  else
    {
      ifstream file(filename.c_str());
      if(!file.good()) throw TigrString("Error opening file: ")+filename;
      loadFromStream(file);
    }
}



void TigrChi2Table::loadFromStream(istream &file)
{
  int row=0;
   
  TigrString line;
  line.getline(file);
  TigrVector<TigrString> &fields=*line.getFields();
  int n=fields.size();
  for(int i=1 ; i<n ; ++i)
    alphas[i-1]=fields[i].asFloat();
  delete &fields;
   
  while(!file.eof())
    {
      line.getline(file);
      if(file.eof()) break;
      TigrVector<TigrString> *fieldsPtr=line.getFields();
      if(!fieldsPtr) continue;
      TigrVector<TigrString> &fields=*fieldsPtr;
      int n=fields.size();
      int df=fields[0].asInt();
      DFs[row]=df;
      for(int i=1 ; i<n ; ++i)
	{
	  float chi=fields[i].asFloat();
	  table[i-1][row]=chi;
	}
      delete fieldsPtr;
      ++row;
    }
}

/*
{
  int row=0;
  TigrFile file(filename);
   
  TigrString line=file.readLine();
  TigrVector<TigrString> &fields=*line.getFields();
  int n=fields.size();
  for(int i=1 ; i<n ; ++i)
    alphas[i-1]=fields[i].asFloat();
  delete &fields;
   
  while(!file.eof())
    {
      TigrString line=file.readLine();
      if(file.eof()) break;
      TigrVector<TigrString> *fieldsPtr=line.getFields();
      if(!fieldsPtr) continue;
      TigrVector<TigrString> &fields=*fieldsPtr;
      int n=fields.size();
      int df=fields[0].asInt();
      DFs[row]=df;
      for(int i=1 ; i<n ; ++i)
	{
	  float chi=fields[i].asFloat();
	  table[i-1][row]=chi;
	}
      delete fieldsPtr;
      ++row;
    }
}
*/


float TigrChi2Table::lookupChi(int df,float alpha)
{
  int alphaIndex=alphaToIndex(alpha);
  int dfIndex=dfToIndex(df);

  if(df<DFs[dfIndex])
    {
      float y1=table[alphaIndex][dfIndex-1];
      float y2=table[alphaIndex][dfIndex];
      int x1=DFs[dfIndex-1];
      int x2=DFs[dfIndex];
      float ratio=float(df-x1)/(x2-x1);
      float chi=y1+ratio*(y2-y1);
      return chi;
    }
  else if(df>DFs[dfIndex] && dfIndex<DFs.size())
    {
      float y1=table[alphaIndex][dfIndex];
      float y2=table[alphaIndex][dfIndex+1];
      int x1=DFs[dfIndex];
      int x2=DFs[dfIndex+1];
      float ratio=float(df-x1)/(x2-x1);
      float chi=y1+ratio*(y2-y1);
      return chi;
    }
  if(dfIndex>=DFs.size()) dfIndex=DFs.size()-1;
  
  return table[alphaIndex][dfIndex];
}



float TigrChi2Table::lookupP(int df,float chi)
{
  int dfIndex=dfToIndex(df);
  unsigned begin=0, end=alphas.size();
  while(begin<end)
    {
      unsigned mid=unsigned((begin+end)/2);
      float midChi=table[mid][dfIndex];
      if(chi>midChi) begin=mid+1;
      else end=mid;
    }
  if(begin>=alphas.size()) begin=alphas.size()-1;
  if(chi<table[begin][dfIndex]) 
    {
      if(begin==0) return alphas[0];
      float x1=table[begin-1][dfIndex];
      float x2=table[begin][dfIndex];
      float y1=alphas[begin-1]; 
      float y2=alphas[begin];
      float ratio=(chi-x1)/float(x2-x1);
      return y1-ratio*(y1-y2);
    }
  else if(chi>table[begin][dfIndex] && begin<alphas.size()-1)
    {
      float x1=table[begin][dfIndex];
      float x2=table[begin+1][dfIndex];
      float y1=alphas[begin];
      float y2=alphas[begin+1];
      float ratio=(chi-x1)/float(x2-x1);
      return y1-ratio*(y1-y2);
    }

  return alphas[begin];
}



unsigned TigrChi2Table::alphaToIndex(float alpha)
{
  unsigned begin=0, end=alphas.size();
  while(begin<end)
    {
      unsigned mid=unsigned((begin+end)/2);
      float midAlpha=alphas[mid];
      if(alpha<midAlpha) begin=mid+1;
      else end=mid;
    }
  if(alphas[begin]!=alpha) 
    throw TigrString("Bad alpha in TigrChi2Table::lookupChi(")+
      alpha+")";
  return begin;
}



unsigned TigrChi2Table::dfToIndex(int df)
{
  unsigned begin=0, end=DFs.size();
  while(begin<end)
    {
      unsigned mid=unsigned((begin+end)/2);
      int midElem=DFs[mid];
      if(df>midElem) begin=mid+1;
      else end=mid;
    }
  return begin;
}
