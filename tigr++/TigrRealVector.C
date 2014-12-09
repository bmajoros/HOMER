#include <math.h>
#include <iostream>
#include "TigrRealVector.H"
using namespace std;



TigrRealVector::TigrRealVector(int dimension)
  : theMatrix(1,dimension)
{
}



double TigrRealVector::operator[](int i) const
{
  return theMatrix(0,i);
}



double &TigrRealVector::operator[](int i)
{
  return theMatrix(0,i);
}



double TigrRealVector::distanceTo(const TigrRealVector &other)
{
  int dim=theMatrix.getNumColumns();
  double sumSquares=0;
  for(int i=0 ; i<dim ; ++i)
    {
      double d=matrix(0,i)-other[i];
      sumSquares+=(d*d);
    }
  return sqrt(sumSquares);
}



double TigrRealVector::dotProduct(const TigrRealVector &other)
{
  double sumProducts=0;
  int columns=theMatrix.getNumColumns();
  for(int i=0 ; i<columns ; ++i)
    sumProducts+=(other.theMatrix(0,i)*theMatrix(0,i));
  return sumProducts;
}



double TigrRealVector::getAngle(const TigrRealVector &other)
{
  double dotProd=dotProduct(other);
  double norm=getMagnitude();
  double normalizedDotProd=dotProd/(other.getMagnitude()*norm);
  double theta=acos(normalizedDotProd);
  return theta;
}



double TigrRealVector::getMagnitude()
{
  double sumSquares=0;

  int cols=theMatrix.getNumColumns();
  for(int i=0 ; i<cols ; ++i)
    {
      double entry=theMatrix(0,i);
      sumSquares+=(entry*entry);
    }

  return sqrt(sumSquares);
}



int TigrRealVector::getDimension()
{
  return theMatrix.getNumColumns();
}



void TigrRealVector::scale(double factor)
{
  int cols=matrix.getNumColumns();
  for(int i=0 ; i<cols ; ++i)
    theMatrix(0,i)*=factor;
}



void TigrRealVector::unitize()
{
  double factor=1.0/getMagnitude();
  scale(factor);
}
