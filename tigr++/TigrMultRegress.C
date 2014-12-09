#include <iostream>
#include "TigrMultRegress.H"
using namespace std;



bool TigrMultRegress::inversionRegression(const TigrDblMatrix &Y,
					  const TigrDblMatrix &X,
					  TigrDblMatrix &coef)
{
  /*
    This method utilizes the fact that (due to matrix algebra):
       XA=Y  ->  A=inv(Xtr*X)*Xtr*Y
    for attribute matrix X, category matrix Y, and coefficient matrix A.
   */

  int numPoints=Y.getNumRows(), numPredictors=X.getNumColumns();
  if(Y.getNumColumns()!=1) 
    throw "TigrMultRegress::regress(): Y matrix must have exactly 1 column";
  if(X.getNumRows()!=numPoints)
    throw "TigrMultRegress::regress(): matrices must have same # of rows";
  TigrDblMatrix Xtr(numPredictors,numPoints);
  X.transpose(Xtr);
  TigrDblMatrix XtrX(numPredictors,numPredictors);
  Xtr.times(X,XtrX);
  TigrDblMatrix XtrXinv(numPredictors,numPredictors);
  if(!XtrX.invert(XtrXinv)) 
    return false;
  TigrDblMatrix XtrXinvXtr(numPredictors,numPoints);
  XtrXinv.times(Xtr,XtrXinvXtr);
  XtrXinvXtr.times(Y,coef);
  return true;
}



bool TigrMultRegress::regress(const TigrDblMatrix &Y,const TigrDblMatrix &X,
			      TigrDblMatrix &coef)
{
  int numColumns=X.getNumColumns(), numRows=X.getNumRows();
  TigrDblMatrix augmentedX(numRows,numColumns+1);
  for(int i=0 ; i<numRows ; ++i)
    {
      for(int j=0 ; j<numColumns ; ++j)
	augmentedX(i,j)=X(i,j);
      augmentedX(i,numColumns)=1.0;
    }
  bool success=inversionRegression(Y,augmentedX,coef);
  return success;
}
