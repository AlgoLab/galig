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
    string real_edges = argv[3];
    string added_edges = argv[4];
    int L = atoi(argv[5]);
    int K = atoi(argv[6]);
    string out_file = argv[7];
    float perc = (100.0-2*L)/100.0;
    cout << "\tStarting..." << endl;
    MemsReader mr = MemsReader(mems);
    mr.readMEMsFile();
    ReferenceGraph g (e_lens, real_edges, added_edges);
    ofstream outFile;
    outFile.open(out_file);

    int i = 0;
    while(mr.hasPattern()) {
	i++;
	pair<string, MemsList> p = mr.popPattern();
	MemsGraph mg (g, p.second, K, perc);
	//mg.saveImage("./out/" + p.first);
	mg.visit();
	mg.saveOutput(outFile, p.first);
    }
    outFile.close();
    cout << "\tEnding." << endl;
}
