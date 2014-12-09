#include <time.h>
#include <sstream>
#include "TigrProgress.H"
using namespace std;



TigrProgress::TigrProgress()
{
  // constructor
}



TigrString TigrProgress::getProgress(unsigned long workDone)
  {
    if(workDone==0) {return TigrString("0% done");}

    float percentDone=
      int(1000*workDone/(float)totalWork+0.50)/10.0;

    unsigned long workRemaining=totalWork-workDone;
    long elapsedSec=time(NULL)-startTime;
    float secPerUnitWork=elapsedSec/(float)workDone;
    int sec=int(secPerUnitWork*workRemaining);

    stringstream os;
    os << percentDone << "% done, ";
    if(sec<60)
      { 
	os << sec << " sec"; 
      }
    else
      {
	float min=int(sec/6.0)/10.0;
	if(min<60)
	  { os << min << " min"; }
	else
	  {
	    float hours=int(min/6.0)/10.0;
	    if(hours<24)
	      { os << hours << " hours"; }
	    else
	      {
		float days=int(hours/2.4)/10.0;
		os << days << " days";
	      }
	  }
      }

    return os.str();
  }



void TigrProgress::start(unsigned long t)
{
  totalWork=t;
  startTime=time(NULL);
}
