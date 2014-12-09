#include <iostream>
#include "TigrInversionMatrix.H"
using namespace std;



TigrInversionMatrix::TigrInversionMatrix(int numRows)
  : TigrDblMatrix(numRows,numRows*2)
{
  setupIdentMatrix();
}



void TigrInversionMatrix::setupIdentMatrix()
{
  TigrDblMatrix &self=*this;
  int rows=getNumRows();
  for(int i=0 ; i<rows ; ++i)
    for(int j=0 ; j<rows ; ++j)
      if(i==j) 
	self.TigrDblMatrix::operator()(i,j)=1.0; 
      else 
	self.TigrDblMatrix::operator()(i,j)=0.0;
}



void TigrInversionMatrix::getInverted(TigrDblMatrix &mat)
{
  TigrDblMatrix &self=*this;
  int rows=getNumRows();
  for(int i=0 ; i<rows ; ++i)
    for(int j=0 ; j<rows ; ++j)
      mat(i,j)=self.TigrDblMatrix::operator()(i,j);
}



double &TigrInversionMatrix::operator()(int row,int column)
{
  int adjustedColumn=TigrDblMatrix::getNumRows()+column;
  return TigrDblMatrix::operator()(row,adjustedColumn);
}



bool TigrInversionMatrix::detectNonInvertible()
{
  TigrDblMatrix &self=*this;
  int rows=getNumRows();
  for(int i=0 ; i<rows ; ++i)
    {
      int j;
      for(j=0 ; j<rows ; ++j)
	if(0.0!=self(i,j)) break;
      if(rows==j) return true;
    }
  return false;
}



int TigrInversionMatrix::getNumColumns()
{
  return TigrDblMatrix::getNumColumns()/2;
}



void TigrInversionMatrix::install(TigrDblMatrix &from)
{
  TigrInversionMatrix &self=*this;
  int rows=getNumRows();
  for(int i=0 ; i<rows ; ++i)
    for(int j=0 ; j<rows ; ++j)
      self(i,j)=from(i,j);
}
