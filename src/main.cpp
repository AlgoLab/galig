/**x
   type* name; -> pointer declaration
   type& name; -> reference declaration

   &var -> pointer to var
   *pointer_var -> var
   **/

#include <iostream>

#include "MEMsReader.hpp"
#include "ReferenceGraph.hpp"
#include "MEMsGraph.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    string mems = argv[1];
    string e_lens = argv[2];
    string edges = argv[3];
    int K = atoi(argv[4]);
    
    cout << "Starting..." << endl;
    MemsReader mr = MemsReader(mems);
    mr.readMEMsFile();
    //mr.print();
    ReferenceGraph g (e_lens, edges);

    while(mr.hasPattern()) {
	pair<string, MemsList> p = mr.popPattern();
	cout << p.first << endl;
	MemsGraph mg (g, p.second, K);
	//mg.save();
	mg.visit();
    }
    cout << "Ending." << endl;
}
