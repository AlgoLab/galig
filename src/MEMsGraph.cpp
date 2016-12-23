#include "MEMsGraph.hpp"

TStr toTStr(const std::string& s) {
  return TStr(s.c_str());
}

std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& vect) {
    for(std::string s : vect) {
	os << s << " ";
    }
    os << "\n";
    return os;
}

TStr MemsGraph::toTStr(const std::string& s) {
  return TStr(s.c_str());
}

int MemsGraph::getNodeId(const std::string& mem) {
    int index;
    try {
	index = MemToIndex.at(mem);
    } catch (const std::out_of_range& oor) {
	index = -1;
    }
    return index;
}

bool MemsGraph::isNode(Mem m) {
    try {
	MemToIndex.at(m.toStr());
	return true;
    }
    catch (const std::out_of_range& oor) {
	return false;
    }
}

void MemsGraph::addNode(Mem mem) {
    std::string label = mem.toStr();
    IndexToMem.insert({nodes_index, mem});
    MemToIndex.insert({label, nodes_index});
    Graph->AddNode(nodes_index);
    labels.AddDat(nodes_index, toTStr(label));
    nodes_index++;
}

void MemsGraph::addEdge(Mem mem1, Mem mem2) {
    Graph->AddEdge(MemToIndex.at(mem1.toStr()), MemToIndex.at(mem2.toStr()));
}

MemsGraph::MemsGraph(ReferenceGraph& g, MemsList& ml, const int& K) {
    Graph = TNGraph::New();
    addNode(Mem(0,0,0));
    int curr_p = 1;
    plen = ml.getLength();
    this->K = K;
    while(curr_p < plen) {
	std::forward_list<Mem> mems1 = ml.getMems(curr_p);
	for(auto it1=mems1.begin(); it1!=mems1.end(); ++it1) {
	    Mem m1 = (*it1);
	    if(!isNode(m1)) {
	     	addNode(m1);
	    	addEdge(Mem(0,0,0), m1);
	    }
	    int i = m1.p + 1;
	    while(i < plen && i < m1.p + m1.l + K) {
		std::forward_list<Mem> mems2 = ml.getMems(i);
		for(auto it2=mems2.begin(); it2!=mems2.end(); ++it2) {
		    Mem m2 = (*it2);
		    //std::cout << "Checking " << m1.toStr() << " -> "<< m2.toStr() << std::endl;
		    if(m1.p + m1.l < m2.p + m2.l) {
			if(m1.t != m2.t && m1.t + m1.l != m2.t + m2.l) {
			    if(g.rank(m1.t - 1) == g.rank(m2.t - 1)) {
				//Stesso esone
				//std::cout << "1" << std::endl;
				if(m2.t > m1.t && m2.t < m1.t + m1.l + K && m1.t + m1.l != m2.t + m2.l) {
				    if(!isNode(m2)) {
				    	//std::cout << "1 Adding " << m2.toStr() << std::endl;
				    	addNode(m2);
				    }
				    /**
				    int wt = m2.t - m1.t - m1.l;
				    int wp = m2.p - m1.p - m1.l;
				    int w;
				    if(wt<0 || wp<0) {
				    	w = abs(wt - wp);
				    }
				    else {
				    	w = max(wt, wp);
				    }
				    **/
				    //std::cout << "1 Adding " << m1.toStr() << " -> " << m2.toStr() << std::endl;
				    addEdge(m1, m2);
				}
			    } else {
				//Esoni diversi
				//std::cout << "2" << std::endl;
				std::vector<int> curr_edge { g.rank(m1.t-1), g.rank(m2.t-1) };
				if(g.contain(curr_edge)) {
				    //std::cout << "." << std::endl;
				    if(m1.t + m1.l >= g.select(g.rank(m1.t-1) + 1) - K && m2.t <= g.select(g.rank(m2.t-1)) + K) {
					if(!isNode(m2)) {
					    //std::cout << "2 Adding " << m2.toStr() << std::endl;
					    addNode(m2);
					}
					/**
					int wt = (g.select(g.rank(m1.t-1) + 1) - m1.t - m1.l) + (m2.t - g.select(g.rank(m2.t-1)) - 1);
					int wp = abs(m2.p - m1.p - m1.l);
					int w = max(wt, wp);
					**/
					addEdge(m1, m2);
					//std::cout << "2 Adding " << m1.toStr() << " -> " << m2.toStr() << std::endl;
				    }
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

    subpaths = std::vector<std::vector<std::vector<int> > >(Graph->GetNodes(), { std::vector<std::vector<int> > { std::vector<int> { } } });
}

void MemsGraph::visit() {
    paths = rec_visit(Graph->BegNI());
}

std::vector<std::vector<int> > MemsGraph::rec_visit(const TNGraph::TNodeI node) {
    int node_id = node.GetId();
    if(subpaths[node_id][0].size() != 0) {
	return subpaths[node_id];
    }
    int out = node.GetOutDeg();
    if(out == 0) {
	std::vector<std::vector<int> > starting_paths { std::vector<int> { node_id } };
	subpaths.at(node_id) = starting_paths;
	return starting_paths;
    }
    int i = 0;
    std::vector<std::vector<int> > paths;
    while(i < out) {
	int child_id = node.GetOutNId(i);
	i++;
	TNGraph::TNodeI child = Graph->GetNI(child_id);
	std::vector<std::vector<int> > starting_paths = rec_visit(child);

	for(std::vector<int> sp : starting_paths) {
	    std::vector<int> p { node_id };
	    p.insert(p.end(),sp.begin(),sp.end());
	    paths.push_back(p);
	}
    }
    subpaths.at(node_id) = paths;
    return paths;
}

void MemsGraph::saveOutput(std::ostream& os, std::string p, const float& perc) {
    for(std::vector<int> path : paths) {
	std::vector<std::string> path_s (path.size(), "");
	Mem starting_mem = IndexToMem.at(path[1]);
	Mem ending_mem = IndexToMem.at(path[path.size()-1]);
	if(ending_mem.p + ending_mem.l - starting_mem.p >= perc*plen) {
	    unsigned int i = 1;
	    while(i<path.size()) {
		path_s[i-1] = labels.GetDat(labels.GetKey(labels.GetKeyId(path[i]))).GetCStr();
		i++;
	    }
	    os << "--- " << p << "\n";
	    os << path_s;
	    os << "\n";
	}
    }
}

void MemsGraph::saveImage(const std::string& patt) {
    std::ofstream myfile;
    myfile.open(patt + ".dot");
    
    std::string dot = "digraph G {\n graph [splines=true overlap=false]\n node  [shape=ellipse, width=0.3, height=0.3]\n";
    for (TNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) { 
	dot += " " + std::to_string(NI.GetId()) + " [label=\"" + labels.GetDat(labels.GetKey(labels.GetKeyId(NI.GetId()))).GetCStr() + "\"];\n";
    }
    
    for (TNGraph::TEdgeI EI = Graph->BegEI(); EI < Graph->EndEI(); EI++) { 
	dot += " " + std::to_string(EI.GetSrcNId()) + " -> " + std::to_string(EI.GetDstNId()) + ";\n";
    }
    dot += "}";
    
    myfile << dot;
    myfile.close();
    if(system(("dot -Tpng ./" + patt + ".dot -o ./" + patt + ".png").c_str()) != 0) {
	std::cerr << "System call error" << std::endl;
    }
}
