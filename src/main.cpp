/**
   type* name; -> pointer declaration
   
   &var -> pointer to var
   *pointer_var -> var
   **/

#include <iostream>

//#include "MEMsList.hpp"
#include "MEMsReader.hpp"

using namespace std;

int main() {
    string p_id = "1";
    int p_len = 5;
    forward_list<vector<int> > mems = {{1,2,3}, {3,2,1}, {3,1,1}, {3,4,1}, {3,5,1}};

    MemsReader mr = MemsReader();
    mr.addPattern(p_id, p_len, mems);
    mr.getPatterns();
    /**
    std::forward_list<std::pair<std::string, MemsList*> > ps = mr.getPatterns();
    for(forward_list<std::pair<std::string, MemsList*> >::iterator it=ps.begin(); it != ps.end(); ++it) {
	std::pair<std::string, MemsList*> p = (*it);
	MemsList ml = *(p.second);
	int i = 0;
	while(i<ml.getLength()) {
	    forward_list<Mem> l = ml.getMems(i);
	    cout << i << endl;
	    for(forward_list<Mem>::iterator it1=l.begin(); it1 != l.end(); ++it1) {
		Mem m = (*it1);
		cout << m.toStr() << endl;
	    }
	    i++;
	}
    }
    **/
/**
    MemsList ml (5);
    ml.addMem(1,2,3);
    ml.addMem(3,2,1);
    ml.addMem(3,1,1);
    ml.addMem(3,4,1);
    ml.addMem(3,5,1);
    int i = 0;
    while(i<ml.getLength()) {
	forward_list<Mem> l = ml.getMems(i);
        if(!l.empty()) {
	    cout << i << endl;
	    for(forward_list<Mem>::iterator it=l.begin(); it != l.end(); ++it) {
		Mem m = (*it);
		cout << m.toStr() << endl;
	    }
	}
	i++;
	}**/
}
