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
    string path = argv[1];
    cout << "Starting..." << endl;
    MemsReader mr = MemsReader(path + "mems");
    mr.readMEMsFile();
    //mr.print();
    ReferenceGraph g (path + "e_lens", path + "edges");

    while(mr.hasPattern()) {
	pair<string, MemsList> p = mr.popPattern();
	cout << p.first << endl;
	MemsGraph mg (g, p.second, 5);
	//mg.save();
	mg.visit();
    }
    cout << "Ending." << endl;
}
