#include <iostream>
#include "TigrStacktrace.H"
using namespace std;



TigrStacktrace::TigrStacktrace(const char *message)
{
  cout << endl << message << endl << endl;
  cout << "Initiating core dump...use your debugger to get a stack trace"
       << endl;
  ((Crasher*)0)->crash();
}



TigrStacktrace::TigrStacktrace(const TigrString &message)
{
  cout << endl << message << endl << endl;
  cout << "Initiating core dump...use your debugger to get a stack trace"
       << endl;
  ((Crasher*)0)->crash();
}
