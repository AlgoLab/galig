#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>


#include <sdsl/bit_vectors.hpp>
#include "Snap.h"

using namespace std;
using namespace sdsl;

typedef struct {
    unsigned int t;
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
	cout << "Unable to open MEMs file" << endl;
    }  
    return MEMs;
}

int getId(TPt<TNodeEDatNet<TInt, TInt> > Graph, TIntStrH labels, mem m) {
    int index = -1;
    for (TNodeEDatNet<TInt, TInt>::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
	if(labels.GetDat(labels.GetKey(labels.GetKeyId(NI.GetId()))).EqI(toTStr(m.toStr()))) {
	    index = NI.GetId();
	    break;
	}
    }
    return index;
}

void saveGraph(TPt<TNodeEDatNet<TInt, TInt> >  Graph, TIntStrH labels) {
    ofstream myfile;
    myfile.open("./out/graph.dot");  
    string dot = "digraph G {\n graph [splines=true overlap=false]\n node  [shape=ellipse, width=0.3, height=0.3]\n";
    for (TNodeEDatNet<TInt, TInt>::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) { 
	dot += " " + to_string(NI.GetId()) + " [label=\"" + labels.GetDat(labels.GetKey(labels.GetKeyId(NI.GetId()))).GetCStr() + "\"];\n";
    }
    for (TNodeEDatNet<TInt, TInt>::TEdgeI EI = Graph->BegEI(); EI < Graph->EndEI(); EI++) { 
	dot += " " + to_string(EI.GetSrcNId()) + " -> " + to_string(EI.GetDstNId()) + " [label=\" " + to_string(EI.GetDat()) + "\"];\n";
    }
    dot += "}";
    myfile << dot;
    myfile.close();
    int r = system("neato -Tpng ./out/graph.dot -o ./out/graph.png");
    cout << "Save graph completed: " << r << endl;
}

vector<int> getExonsLengths(string fpath) {
    vector<int> e_lens;
    string line;
    ifstream memsFile(fpath);
    if (memsFile.is_open()) {
	while(getline(memsFile,line)) {
	    e_lens.push_back(stoi(line));
	}
	memsFile.close();
    }
    else {
	cout << "Unable to open exons file" << endl;
    }  
    return e_lens;
}

vector<vector<int> > DFS(TPt<TNodeEDatNet<TInt, TInt> > Graph, TIntStrH labels, TNodeEDatNet<TInt, TInt>::TNodeI node) {
    vector<vector<int> > paths;  
    /**
       int out = node.GetOutDeg();
       int node_id = node.GetId();
       int i = 0;
       cout << "\t--------------\nNode: " << node.GetId() << endl;
       cout << "Out: " << out << endl;
       cout << labels.GetDat(labels.GetKey(labels.GetKeyId(node.GetId()))).GetCStr() << endl;
       while(i < out) {
       int child_id = node.GetOutNId(i);
       cout << "Child: " << child_id << endl;
       cout << labels.GetDat(labels.GetKey(labels.GetKeyId(child_id))).GetCStr() << endl;
       TNodeEDatNet<TInt, TInt>::TNodeI child = Graph->GetNI(child_id);
       TNodeEDatNet<TInt, TInt>::TEdgeI edge = Graph->GetEI(node_id, child_id);
       cout << to_string(edge.GetDat()) << endl;
       vector<vector<int> > subpaths = DFS(Graph, labels, child);
       for(auto sp : subpaths) {
       vector<int> p { node_id };
       p.insert(p.end(),sp.begin(),sp.end());
       for(auto n : p) {
       cout << "\t\t" << n << " ";
       }
       paths.push_back(p);
       }
       i++;
       }
  
       for(auto p : paths) {
       for(auto n : p) {
       cout << n << " ";
       }
       cout << endl << endl;
       }
    **/
    return paths;
}


int main(int argc, char* argv[]) {
    // Input
    string mems_file = argv[1];
    string e_lens_file = argv[2];
  
    vector<int> e_lens = getExonsLengths(e_lens_file);
    int L = stoi(argv[3]);
    int k = stoi(argv[4]);
    //text_path
    //pattern_path
  
    int plen = stoi(argv[5]);
    //int e_lens[] {4,5};
  
    int tot_L = 1;
    for(int l:e_lens) {
	tot_L += l+1;
    }
  
    //Bit Vector Setup
    bit_vector BV(tot_L, 0);
  
    int i = 0;
    BV[i] = 1;
    for(int l:e_lens) {
	i += l+1;
	BV[i] = 1;
    }
    i = 0;
    while(i<tot_L) {
	cout << BV[i];
	i++;
    }
    cout << endl;
    /**
     * Template: <uint16_t t_bs = 63, class t_rac = int_vector<>, uint16_t t_k = 32>
     *      - t_bs: Size of a basic block. (127)
     *      - t_rac: Random access integer vector. Use to store the block types.
     *      - t_k: A rank sample value is stored before every t_k-th basic block.
     **/
    rrr_vector<> rrrb(BV);
    rrr_vector<>::rank_1_type rank_BV(&rrrb); //rank_BV(int i)
    rrr_vector<>::select_1_type select_BV(&rrrb); //select_BV(int i)
  
    //Extracting MEMs from file
    vector<vector<mem > > MEMs = extractMEMs(mems_file, plen);
  
    //Printing MEMs
    int p_ = 1;
    while(p_ < plen) {
	for(mem m : MEMs[p_]) {
	    cout << m.toStr() << endl;
	}
	p_++;
    }
  
    //MEMsGraph
    TPt<TNodeEDatNet<TInt, TInt> >  Graph = TNodeEDatNet<TInt, TInt>::New();
    TIntStrH labels;
    int nodes_index = 1;
    //int curr_index;
  
    Graph->AddNode(0, 0);
    labels.AddDat(0, "Start");
  
    unsigned int curr_p = 1;
    while(curr_p<MEMs.size()) {
	for(mem m1 : MEMs[curr_p]) {
	    int m1_index = getId(Graph, labels, m1);
      
	    if(m1_index == -1) {
		m1_index = nodes_index;
		nodes_index++;
		Graph->AddNode(m1_index, m1.l);
		Graph->AddEdge(0, m1_index);
		labels.AddDat(m1_index, toTStr(m1.toStr()));
	    }
      
	    int i = m1.p + 1;
	    while(i < plen and i < m1.p + m1.l + k) {
		for(mem m2 : MEMs[i]) { //Per tutti i m2 "consecutivi" a m1
		    if(m1.p + m1.l != m2.p + m2.l) { //Se m1 e m2 non finiscono nello stesso punto sul pattern
			if(m1.t != m2.t && m1.t + m1.l != m2.t + m2.l) { //Se m1 e m2 non iniziano e finiscono negli stessi punti sul testo
			    if(rank_BV(m1.t - 1) == rank_BV(m2.t - 1)) { //Se m1 e m2 sono nello stesso nodo
				cout << "(1) Checking " << m1.toStr() << " -> " << m2.toStr() << endl;
				if(m2.t > m1.t && m2.t < m1.t + m1.l + k && m1.t + m1.l != m2.t + m2.l) {
				    cout << "\tLinking " << m1.toStr() << " to " << m2.toStr() << endl;
				    int m2_index = getId(Graph, labels, m2);
                  
				    if(m2_index == -1) {
					m2_index = nodes_index;
					nodes_index++;
					Graph->AddNode(m2_index, m2.l);
					labels.AddDat(m2_index, toTStr(m2.toStr()));
				    }
                  
				    //Weight
				    int wt = m2.t - m1.t - m1.l;
				    if(wt<0) {
					wt = abs((m1.t + m1.l - m2.t) - (m1.p + m1.l - m2.p));
				    }
                  
				    int wp = m2.p - m1.p - m1.l;
				    if(wp<0) {
					wp = abs((m1.t + m1.l - m2.t) - (m1.p + m1.l - m2.p));
				    }
                  
				    int w = max(wt, wp);
                  
				    Graph->AddEdge(m1_index, m2_index, w);
				}
			    }
			    else { //Se m1 e m2 sono in due nodi differenti
				cout << "(2) Checking " << m1.toStr() << " -> " << m2.toStr() << endl;
				/**
				 * for debug
				 int x1 = rank_BV(m1.t-1);
				 int y1 = select_BV(x1 + 1);
				 int x2 = rank_BV(m2.t-1);
				 int y2 = select_BV(x2 + 1);
				*/
				if(m1.t + m1.l >= select_BV(rank_BV(m1.t-1) + 1) - k && m2.t <= select_BV(rank_BV(m2.t-1)) + k) {
				    cout << "\tLinking " << m1.toStr() << " to " << m2.toStr() << endl;
				    int m2_index = getId(Graph, labels, m2);
                
				    if(m2_index == -1) {
					m2_index = nodes_index;
					nodes_index++;
					Graph->AddNode(m2_index);
					labels.AddDat(m2_index, toTStr(m2.toStr()));
				    }
				    //Weight
				    int wt = (select_BV(rank_BV(m1.t-1) + 1) - m1.t - m1.l) + (m2.t - select_BV(rank_BV(m2.t-1)) - 1);
                  
				    int wp = abs(m2.p - m1.p - m1.l);
                  
				    int w = max(wt, wp);
				    Graph->AddEdge(m1_index, m2_index, w);
				}
			    }
			}
		    }
		}
		i++;
	    }
	}
	curr_p++;
    }
    Graph->AddNode(nodes_index, 0);
    labels.AddDat(nodes_index, "End");
    for (TNodeEDatNet<TInt, TInt>::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) { 
	if(NI.GetOutDeg() == 0 && NI.GetId() != nodes_index) {
	    Graph->AddEdge(NI.GetId(), nodes_index, 0);
	}
    }
  
    saveGraph(Graph, labels);
  
    //Visits
    vector<vector<int> > paths = DFS(Graph, labels, Graph->BegNI());
  
    cout << paths.size() << endl;
  
    for(auto p : paths) {
	for(auto n : p) {
	    cout << n << " ";
	}
	cout << endl << endl;
    }
  
    //PNGraph tree = TSnap::GetBfsTree(Graph, 0, true, true);
    //TSnap::DrawGViz(tree, gvlDot, "BFS.png", "");
}
