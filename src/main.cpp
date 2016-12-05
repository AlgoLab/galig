#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "Snap.h"


using namespace std;


const string IO_folder = "../IO/";

typedef struct {
  int t;
  int p;
  int l;
  
  string toStr() {
    return "(" + to_string(t) + "," + to_string(p) + "," + to_string(l) + ")";
  }
} mem;

TStr toTStr(string s) {
  return TStr(s.c_str());
}

mem extractMEM(string line) {
  mem m;
  string delimiter = ",";
  size_t pos = 0;
  string token;
  bool flag = true;
  while ((pos = line.find(delimiter)) != string::npos) {
    token = line.substr(0, pos);
    if(flag) {
      m.t = stoi(token);
    }
    else {
      m.p = stoi(token);
    }
    line.erase(0, pos + delimiter.length());
    flag = not(flag);
  }
  m.l = stoi(line);

  return m;
}

vector<vector<mem > > extractMEMs(string fpath, int plen) {
  vector<vector<mem > > MEMs (plen + 1, vector<mem >());

  string line;
  ifstream memsFile(fpath);
  if (memsFile.is_open()) {
    while(getline(memsFile,line)) {
      mem m = extractMEM(line);
      MEMs[m.p].push_back(m);
    }
    memsFile.close();
  }
  else {
    cout << "Unable to open file" << endl;
  }
  
  return MEMs;
}

int getId(PNGraph Graph, TIntStrH labels, mem m) {
  int index = -1;
  for (TNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    if(labels.GetDat(labels.GetKey(labels.GetKeyId(NI.GetId()))).EqI(toTStr(m.toStr()))) {
      index = NI.GetId();
      break;
    }
  }
  return index;
}

int main() {
  int L = 2;
  int plen = 9;
  int k = 3;
  string fpath = IO_folder + "mems";
  vector<vector<mem > > MEMs = extractMEMs(fpath, plen);
  
  int p_ = 1;
  while(p_ < plen) {
    for(mem m : MEMs[p_]) {
      cout << m.toStr() << endl;
    }
    p_++;
  }
  
  PNGraph Graph = TNGraph::New();
  TIntStrH labels;
  int nodes_index = 0;
  int curr_index;
  
  int curr_p = 1;
  while(curr_p<MEMs.size()) {
    for(mem m1 : MEMs[curr_p]) {
      int m1_index = getId(Graph, labels, m1);
      
      if(m1_index == -1) {
        m1_index = nodes_index;
        nodes_index++;
        Graph->AddNode(m1_index);
        labels.AddDat(m1_index, toTStr(m1.toStr()));
      }
      
      int i = m1.p + 1;
      while(i < plen and i < m1.p + m1.l + k) {
        for(mem m2 : MEMs[i]) {
          if(m1.p + m1.l != m2.p + m2.l) {
          //if(m1.id == m2.id) {
            cout << "Checking " << m1.toStr() << " -> " << m2.toStr() << endl;
            if(m2.t > m1.t && m2.t < m1.t + m1.l + k && m1.t + m1.l != m2.t + m2.l) {
              cout << "\tLinking " << m1.toStr() << " to " << m2.toStr() << endl;
              int m2_index = getId(Graph, labels, m2);
              
              if(m2_index == -1) {
                m2_index = nodes_index;
                nodes_index++;
                Graph->AddNode(m2_index);
                labels.AddDat(m2_index, toTStr(m2.toStr()));
              }
              Graph->AddEdge(m1_index, m2_index);
              
            }
          }
          //else {
          // 
          //}
        }
        i++;
      }
    }
    curr_p++;
  }
  //string png_path = IO_folder + "graph.png";
  TSnap::DrawGViz(Graph, gvlDot, toTStr(IO_folder + "graph.png"), "", labels);
}
