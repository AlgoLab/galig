/**
   type* name; -> pointer declaration
   type& name; -> reference declaration

   &var -> pointer to var
   *pointer_var -> var
**/

#include <iostream>
#include <fstream>

#include "MEMsReader.hpp"
#include "ReferenceGraph.hpp"
#include "MEMsGraph.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    string mems = argv[1];
    string e_lens = argv[2];
    string edges = argv[3];
    int K = atoi(argv[4]);
    string out_file = argv[5];

    cout << "Starting..." << endl;
    MemsReader mr = MemsReader(mems);
    mr.readMEMsFile();
    //mr.print();
    ReferenceGraph g (e_lens, edges);

    ofstream outFile;
    outFile.open(out_file);

    while(mr.hasPattern()) {
	pair<string, MemsList> p = mr.popPattern();
	outFile << p.first << "\n";
	MemsGraph mg (g, p.second, K);
	//mg.saveImage("./out/" + p.first);
	mg.visit();
	mg.saveOutput(outFile);
	outFile << "\n-------------------------------------------------------\n";
    }
    outFile.close();
    cout << "Ending." << endl;
}
