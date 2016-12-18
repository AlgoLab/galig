/**
   type* name; -> pointer declaration
   type& name; -> reference declaration

   &var -> pointer to var
   *pointer_var -> var
   **/

#include <iostream>

#include "MEMsReader.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    if(argc == 1) {
	MemsReader mr = MemsReader("./example/tmp/mems");
	mr.readMEMsFile(70);
	//cout << "end" << endl;
	mr.print();
    }
    /**
    if(argc == 3) {
	cout << "Starting..." << endl;
	MemsReader mr = MemsReader(argv[1]);
	mr.readMEMsFile(stoi(argv[2]));
	mr.print();
    }
    **/
}
