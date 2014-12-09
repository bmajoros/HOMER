#include <iostream>
#include "TigrContingencyTbl.H"
using namespace std;



TigrContingencyTbl::TigrContingencyTbl(int width,int height)
  : TigrArray2D<T>(width,height), 
    rowTotals(height), columnTotals(width)
{
  // ctor
}



int TigrContingencyTbl::getColumnTotal(int column)
{
  return columnTotals[column];
}



int TigrContingencyTbl::getGrandTotal()
{
  return grandTotal;
}



int TigrContingencyTbl::getRowTotal(int row)
{
  return rowTotals[row];
}



void TigrContingencyTbl::computeTotals()
{
  const int width=getFirstDimension();
  const int height=getSecondDimension();

  columnTotals.setAllTo(0);
  rowTotals.setAllTo(0);
  grandTotal=0;

  for(int x=0 ; x<width ; ++x)
    for(int y=0 ; y<height ; ++y)
      {
	int entry=(*this)[x][y];
	grandTotal+=entry;
	columnTotals[x]+=entry;
	rowTotals[y]+=entry;
      }
}
