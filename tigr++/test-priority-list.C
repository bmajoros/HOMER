
// test-priority-list.C

#include <string>
#include <iostream>
#include "tigr++/TigrCommandLine.H"
#include "tigr++/FixedSizePriorityList.H"

class Application
{
public:
  Application();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
  {
    try
      {
	Application app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



Application::Application()
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    TigrCommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=0) 
      throw string("test-priority-list <no arguments>");

    FixedSizePriorityList<int> Q(5);
    Q.insert(9);
    Q.insert(6);
    Q.insert(8);
    Q.insert(7);
    Q.insert(5);
    Q.insert(4);
    Q.insert(1);
    Q.insert(2);
    Q.insert(3);
    Q.insert(10);

    Q.erase(Q.begin());

    FixedSizePriorityList<int>::iterator cur=Q.begin(), end=Q.end();
    for(; cur!=end ; ++cur)
      cout<<"visiting "<<*cur<<endl;

    return 0;
  }

