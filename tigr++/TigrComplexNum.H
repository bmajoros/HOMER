/**************************************************************
TigrComplexNum.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrComplexNum_H
#define INCL_TigrComplexNum_H
using namespace std;

class TigrComplexNum
{
public:
  TigrComplexNum(double real=0.0,double imag=0.0);
  TigrComplexNum operator+(const TigrComplexNum &);
  TigrComplexNum operator*(const TigrComplexNum &);
  TigrComplexNum getConjugate();
  double &getImag();
  double &getReal();
  double getImag() const;
  double getModulus() const;
  double getReal() const;
private:
  double real, imag;
};

#endif
