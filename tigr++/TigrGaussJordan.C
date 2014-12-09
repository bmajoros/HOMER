#include <iostream>
#include "TigrGaussJordan.H"
using namespace std;



bool TigrGaussJordan::invert(TigrDblMatrix &thisMatrix,TigrDblMatrix &intoHere)
{
  int numRows=thisMatrix.getNumRows();
  TigrInversionMatrix inversionMatrix(numRows);
  inversionMatrix.install(thisMatrix);

  for(int i=0 ; i<numRows ; ++i)
    {
      int nonzeroRow=getNonzeroRow(i,inversionMatrix);
      if(nonzeroRow<0) return false;
      if(nonzeroRow!=i) 
	inversionMatrix.swapRows(nonzeroRow,i);
      double firstNonzeroEntry=inversionMatrix(i,i);
      if(firstNonzeroEntry!=1.0)
	{
	  double factor=1.0/firstNonzeroEntry;
	  inversionMatrix.multiplyRowBy(i,factor);
	}
      zeroOut(i,inversionMatrix);
    }
  inversionMatrix.getInverted(intoHere);
  return true;
}



int TigrGaussJordan::getNonzeroRow(int startRow,TigrInversionMatrix &inMatrix)
{
  int numRows=inMatrix.getNumRows();
  for(int i=startRow ; i<numRows ; ++i)
    if(inMatrix(i,startRow)!=0.0)
      return i;
  return -1;
}



void TigrGaussJordan::zeroOut(int exceptRow,TigrInversionMatrix &inMatrix)
{
  int numRows=inMatrix.getNumRows(); 
  for(int i=0 ; i<numRows ; ++i)
    if(i!=exceptRow)
      inMatrix.addRowMultiple(exceptRow,i,-inMatrix(i,exceptRow));
}
