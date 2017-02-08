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
  }
  outFile.close();
}
