#include "TigrTime.H"
using namespace std;



TigrString getDateAndTime()
{
  time_t t=time(NULL);
  return TigrString(ctime(&t));
}



TigrString TigrTime::elapsedTime()
{
  // seconds
  float sec=elapsedSeconds();
  if(sec<60)
    return TigrString("")+sec+" sec";

  // minutes
  float min=sec/60;
  if(min<60)
    return TigrString("")+min+" min";

  // hours
  float hours=min/60;
  if(hours<24)
    return TigrString("")+hours+" hours";
  
  // days
  float days=hours/24;
  return TigrString("")+days+" days";
}



float TigrTime::elapsedSeconds()
{
  float seconds=timevalStop.tv_sec-timevalStart.tv_sec;
  float uSec=timevalStop.tv_usec-timevalStart.tv_usec;
  
  return seconds+uSec/1000000.0;
}



void TigrTime::startCounting()
{
  gettimeofday(&timevalStart,NULL);
}



void TigrTime::stopCounting()
{
  gettimeofday(&timevalStop,NULL);
}
