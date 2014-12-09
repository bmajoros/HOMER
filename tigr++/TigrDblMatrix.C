#include <iostream>
#include "TigrStacktrace.H"
#include "TigrGaussJordan.H"
#include "TigrDblMatrix.H"
using namespace std;



TigrDblMatrix::TigrDblMatrix(int rows,int columns) 
  : theArray(rows,columns)
{
}



TigrDblMatrix::TigrDblMatrix(const TigrDblMatrix &m)
  : theArray(m.theArray)
{
}



ostream &operator<<(ostream &os,TigrDblMatrix &mat)
{
  mat.printOn(os);
  return os;
}



double TigrDblMatrix::operator()(int row,int col) const
{
  TigrArray2D<double> &array=const_cast<TigrArray2D<double>&>(theArray);
  return array[row][col];
}



double &TigrDblMatrix::operator()(int row,int col)
{
  TigrArray2D<double> &array=const_cast<TigrArray2D<double>&>(theArray);
  return array[row][col];
}



TigrDblMatrix &TigrDblMatrix::operator=(TigrDblMatrix &other)
{
  int rows=theArray.getFirstDim();
  int cols=theArray.getSecondDim();
  for(int i=0 ; i<rows ; ++i)
    for(int j=0 ; j<cols ; ++j)
      theArray[i][j]=other.theArray[i][j];
}



bool TigrDblMatrix::invert(TigrDblMatrix &resultingMatrix) const
{
  int rows=getNumRows(), cols=getNumColumns();
  if(cols!=rows)
    throw TigrStacktrace("matrix is non-square in TigrDblMatrix::invert()");
  if(resultingMatrix.getNumColumns()!=cols ||
     resultingMatrix.getNumRows()!=rows)
    throw TigrStacktrace("resulting matrix of wrong size in TigrDblMatrix::invert()");

  TigrDblMatrix &self=const_cast<TigrDblMatrix&>(*this);
  bool success=TigrGaussJordan::invert(self,resultingMatrix);
  return success;
}



int TigrDblMatrix::getNumColumns() const
{
  TigrArray2D<double> &array=const_cast<TigrArray2D<double>&>(theArray);
  return array.getSecondDim();
}



int TigrDblMatrix::getNumRows() const
{
  TigrArray2D<double> &array=const_cast<TigrArray2D<double>&>(theArray);
  return array.getFirstDim();
}



void TigrDblMatrix::addRowMultiple(int sourceRow,int destinationRow,
				   double factor)
{
  int cols=theArray.getSecondDim();
  for(int i=0 ; i<cols ; ++i)
    theArray[destinationRow][i]+=theArray[sourceRow][i]*factor;
}



void TigrDblMatrix::getColumn(int col,TigrDblArray1D &into) const
{
  TigrArray2D<double> &array=const_cast<TigrArray2D<double>&>(theArray);
  int rows=getNumRows();
  if(into.size()!=rows)
    throw TigrStacktrace("wrong sized vector in TigrDblMatrix::getColumn()");
  for(int i=0 ; i<rows ; ++i) into[i]=array[i][col];
}



void TigrDblMatrix::getRow(int row,TigrDblArray1D &into) const
{
  TigrArray2D<double> &array=const_cast<TigrArray2D<double>&>(theArray);
  int cols=getNumColumns();
  if(into.size()!=cols)
    throw TigrStacktrace("wrong sized vector in TigrDblMatrix::getRow()");
  for(int i=0 ; i<cols ; ++i)
    into[i]=array[row][i];
}



void TigrDblMatrix::multiplyRowBy(int whichRow,double factor)
{
  int cols=theArray.getSecondDim();
  for(int i=0 ; i<cols ; ++i) theArray[whichRow][i]*=factor;
}



void TigrDblMatrix::printOn(ostream &os)
{
  int rows=theArray.getFirstDim();
  int cols=theArray.getSecondDim();
  for(int i=0 ; i<rows ; ++i)
    {
      for(int j=0 ; j<cols ; ++j)
	os << theArray[i][j] << '\t';
      os << endl;
    }
}



void TigrDblMatrix::setAllTo(double to)
{
  theArray.setAllTo(to);
}



void TigrDblMatrix::swapRows(int r1,int r2)
{
  double tempCell;
  int cols=theArray.getSecondDim();
  for(int i=0 ; i<cols ; ++i)
    {
      tempCell=theArray[r2][i];
      theArray[r2][i]=theArray[r1][i];
      theArray[r1][i]=tempCell;
    }
}



void TigrDblMatrix::times(const TigrDblMatrix &otherMatrix,
			  TigrDblMatrix &resultMatrix) const
{
  TigrDblMatrix &thisOne=*const_cast<TigrDblMatrix*>(this);
  TigrDblMatrix &thatOne=const_cast<TigrDblMatrix&>(otherMatrix);
  multiply(thisOne,thatOne,resultMatrix);
}



void TigrDblMatrix::transpose(TigrDblMatrix &result) const
{
  int nRows=getNumRows(), nCols=getNumColumns();
  if(nCols!=result.getNumRows() ||
     nRows!=result.getNumColumns())
    throw TigrStacktrace("wrong-sized matrix in TigrDblMatrix::transpose");
  TigrArray2D<double> &thatArray=result.theArray;
  TigrArray2D<double> &thisArray=const_cast<TigrArray2D<double>&>(theArray);
  for(int i=0 ; i<nRows ; ++i)
    for(int j=0 ; j<nCols ; ++j)
      thatArray[j][i]=thisArray[i][j];
}



void multiply(TigrDblMatrix &leftM,TigrDblMatrix &rightM,
	      TigrDblMatrix &resultMatrix)
{
  int rCols=resultMatrix.getNumColumns();
  int rRows=resultMatrix.getNumRows();
  int rowLength=leftM.getNumColumns();

  if(resultMatrix.getNumRows()!=leftM.getNumRows() ||
     resultMatrix.getNumColumns()!=rightM.getNumColumns())
    throw TigrStacktrace("wrong-sized matrix in TigrDblMatrix.C");
  if(rightM.getNumRows()!=leftM.getNumColumns())
    throw TigrStacktrace("multiplying differently-sized matrices");

  double s;
  for(int i=0 ; i<rRows ; ++i)
    for(int j=0 ; j<rCols ; ++j)
      {
	s=0.0;
	for(int k=0 ; k<rowLength ; ++k)
	  {
	    double x=rightM(k,j)*leftM(i,k);
	    s+=x;
	  }
	resultMatrix(i,j)=s;
      }
}
