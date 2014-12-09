
// test-binary-tree.C

#include <string>
#include <iostream>
#include "TigrCommandLine.H"
#include "TigrFastBinTree.H"
#include "FixedSizePriorityQueue.H"


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


#include <stl>

int Application::main(int argc,char *argv[])
  {
    // Process command line
    TigrCommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=0) 
      throw string("test-binary-tree <no arguments>");

    class Comp : public TigrComparator<int>
    {
    public:
      virtual bool equal(int &a,int &b) {return a==b;}
      virtual bool greater(int &a,int &b) {return a<b;}
      virtual bool less(int &a,int &b) {return a>b;}
    };

    Comp comp;
    FixedSizePriorityQueue<int> t(5,comp);
    t.insert(6);
    t.insert(9);
    t.insert(3);
    t.insert(1);
    t.insert(7);
    t.insert(10);
    t.insert(8);
    t.insert(5);
    t.insert(4);
    t.insert(2);

    if(!t.isEmpty())
      cout<<"max="<<t.peekMax()<<", min="<<t.peekMin()<<endl;

    TigrFastBinTree<int>::iterator cur=t.begin(), end=t.end();
    for(; cur!=end ; ++cur)
      cout<<"visiting "<<*cur<<endl;
    cout<<endl;

    return 0;
  }

