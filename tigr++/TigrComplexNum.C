#include <iostream>
#include "TigrComplexNum.H"
using namespace std;



TigrComplexNum::TigrComplexNum(double real=0.0,double imag=0.0)
  : real(real), imag(imag)
{
   
}



TigrComplexNum TigrComplexNum::operator+(const TigrComplexNum &other)
{
  return TigrComplexNum(real+other.real,imag+other.imag);
}



TigrComplexNum TigrComplexNum::operator*(const TigrComplexNum &other)
{
  double a=real, b=imag, c=other.real, d=other.imag;
  return TigrComplexNum(a*c-b*d,a*d+c*b);
}



TigrComplexNum TigrComplexNum::getConjugate()
{
  return TigrComplexNum(real,-imag);
}



double &TigrComplexNum::getImag()
{
  return imag;
}



double &TigrComplexNum::getReal()
{
  return real;
}



double TigrComplexNum::getImag() const
{
  return imag;
}



double TigrComplexNum::getModulus() const
{
  return sqrt(real*real + imag*imag);
}



double TigrComplexNum::getReal() const
{
  return real;
}
