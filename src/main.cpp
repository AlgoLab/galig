/**x
   type* name; -> pointer declaration
   type& name; -> reference declaration

   &var -> pointer to var
   *pointer_var -> var
   **/

#include <iostream>

#include "MEMsReader.hpp"
#include "ReferenceGraph.hpp"

using namespace std;

//int main(int argc, char* argv[]) {
int main() {
    cout << "Starting..." << endl;
    MemsReader mr = MemsReader("./example/in/mems");
    mr.readMEMsFile();
    //mr.print();
    ReferenceGraph g ("./example/in/e_lens", "./example/in/edges");
    /**
    if(mr.hasPattern()) {
	pair<string, MemsList> p = mr.popPattern();
	cout << p.first << endl;
    }
    **/
    cout << "Ending." << endl;
}
