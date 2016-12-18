//=================================
// include guard
#ifndef MEM
#define MEM

//=================================
// forward declared dependencies
//class Foo;

//=================================
// included dependencies
#include <string>

class Mem {
private:
    int t;
    int p;
    int l;
public:
    Mem();
    Mem(int t, int p, int l);
    int getT();
    int getP();
    int getL();
    std::string toStr();
};

#endif
