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
    float perc = (100-2*L);

    MemsReader mr = MemsReader(mems);
    mr.readMEMsFile();

    ReferenceGraph g (e_lens, real_edges, added_edges);
    ofstream outFile;
    outFile.open(out_file);

    while(mr.hasPattern()) {
	pair<string, MemsList> p1 = mr.popPattern();
	pair<string, MemsList> p2 = mr.popPattern();

	MemsGraph mg1 (g, p1.second, K, perc);
	MemsGraph mg2 (g, p2.second, K, perc);

	mg1.visit();
	mg2.visit();

        list<string> out1 (mg1.getOutput());
	list<string> out2 (mg2.getOutput());

	if(out1.size() > 0) {
	    for(string s : out1) {
		outFile << p1.first << " " << s << "\n";
	    }
	}
	else if(out2.size() > 0) {
	    for(string s : out2) {
		outFile << p2.first << " " << s << "\n";
	    }
	}
	/**
        if(out1.first >= 0 && out2.first >= 0) {
	    cout << "1" << endl;
	    if(out1.first <= out2.first) {
		outFile << p1.first << " " << out1.second << " " << out1.first << "\n";
	    }
	    else {
		outFile << p2.first << " " << out2.second << " " << out2.first << "\n";
	    }
	}
	else if(out1.first == -1 && out2.first >= 0) {
	    cout << "2" << endl;
	    outFile << p2.first << " " << out2.second << " " << out2.first << "\n";
	}
	else if(out1.first >= 0 && out2.first == -1) {
	    cout << "2" << endl;
	    outFile << p1.first << " " << out1.second << " " << out1.first << "\n";
	}
	**/
    }
    outFile.close();
}
