#include "MEMsGraph.hpp"

std::ostream& operator<<(std::ostream& os, const std::vector<int>& vect) {
    for(int i : vect) {
	os << i << " ";
    }
    os << "\n";
    return os;
}

TStr MemsGraph::toTStr(const std::string& s) {
  return TStr(s.c_str());
}

int MemsGraph::getNodeId(const std::string& mem) {
    int node_index = -1;
    for (TNodeEDatNet<TInt, TInt>::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
	if(labels.GetDat(labels.GetKey(labels.GetKeyId(NI.GetId()))).EqI(toTStr(mem))) {
	    node_index = NI.GetId();
	    break;
	}
    }
    return node_index;
}

void MemsGraph::addNode(const int& exon_index, const std::string& label) {
    Graph->AddNode(nodes_index, exon_index);
    labels.AddDat(nodes_index, toTStr(label));
    nodes_index++;
}

void MemsGraph::addEdge(const int& exon_index_1, const int& exon_index_2, const int& w) {
    Graph->AddEdge(exon_index_1, exon_index_2, w);
}

MemsGraph::MemsGraph(ReferenceGraph &g, MemsList ml, const int& K) {
    Graph = TNodeEDatNet<TInt, TInt>::New();
    addNode(0, "Start");
    
    int curr_p = 1;
    int plen = ml.getLength();
    while(curr_p < plen) {
        std::forward_list<Mem> mems1 = ml.getMems(curr_p);
	for(auto it1=mems1.begin(); it1!=mems1.end(); ++it1) {
	    Mem m1 = (*it1);
	    int m1_index = getNodeId(m1.toStr());
	    if(m1_index == -1) {
		m1_index = nodes_index;
		addNode(m1_index, m1.toStr());
		addEdge(0, m1_index, 0);
	    }

	    int i = m1.p + 1;
	    while(i < plen && i < m1.p + m1.l + K) {
		std::forward_list<Mem> mems2 = ml.getMems(i);
		for(auto it2=mems2.begin(); it2!=mems2.end(); ++it2) {
		    Mem m2 = (*it2);
		    if(m1.p + m1.l != m2.p + m2.l) {
			if(m1.t != m2.t && m1.t + m1.l != m2.t + m2.l) {
			    //std::cout << m1.t-1 << " " << m2.t-1 << std::endl;
			    //std::cout << g.rank(m1.t - 1) << " " << g.rank(m2.t-1) << std::endl;
			    if(g.rank(m1.t - 1) == g.rank(m2.t - 1)) {
				if(m2.t > m1.t && m2.t < m1.t + m1.l + K && m1.t + m1.l != m2.t + m2.l) {
				    int m2_index = getNodeId(m2.toStr());
				    if(m2_index == -1) {
					//std::cout << "1 Adding " << m2.toStr() << std::endl;
					m2_index = nodes_index;
					addNode(m2_index, m2.toStr());
				    }
				    
				    int wt = m2.t - m1.t - m1.l;
				    int wp = m2.p - m1.p - m1.l;
				    int w;
				    if(wt<0 || wp<0) {
					w = abs(wt - wp);
				    }
				    else {
					w = max(wt, wp);
				    }
				    
				    //std::cout << "1 Adding " << m1.toStr() << " -> " << m2.toStr() << std::endl;
				    addEdge(m1_index, m2_index, w);
				}
			    }
			    else {
				std::vector<int> curr_edge { g.rank(m1.t-1), g.rank(m2.t-1) };
				if(g.contain(curr_edge)) {
				    if(m1.t + m1.l >= g.select(g.rank(m1.t-1) + 1) - K && m2.t <= g.select(g.rank(m2.t-1)) + K) {
					int m2_index = getNodeId(m2.toStr());
					if(m2_index == -1) {
					    //std::cout << "2 Adding " << m2.toStr() << std::endl;
					    m2_index = nodes_index;
					    addNode(m2_index, m2.toStr());
					}
					int wt = (g.select(g.rank(m1.t-1) + 1) - m1.t - m1.l) + (m2.t - g.select(g.rank(m2.t-1)) - 1);
					int wp = abs(m2.p - m1.p - m1.l);
					int w = max(wt, wp);
				        addEdge(m1_index, m2_index, w);
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

    int end_index = nodes_index;
    addNode(0, "End");
    for(TNodeEDatNet<TInt, TInt>::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
	if(NI.GetOutDeg() == 0 && NI.GetId() != end_index) {
	    Graph->AddEdge(NI.GetId(), end_index, 0);
	}
    }

    subpaths = std::vector<std::vector<std::vector<int> > >(Graph->GetNodes(), { std::vector<std::vector<int> > { std::vector<int> { } } });
}

void MemsGraph::visit() {
    paths = rec_visit(Graph->BegNI());
}

std::vector<std::vector<int> > MemsGraph::rec_visit(const TNodeEDatNet<TInt, TInt>::TNodeI node) {
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
	TNodeEDatNet<TInt, TInt>::TNodeI child = Graph->GetNI(child_id);
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

void MemsGraph::saveOutput(std::ostream& os) {
    for(std::vector<int> path : paths) {
	os << path;
    }
}

void MemsGraph::saveImage(const std::string& patt) {
    std::ofstream myfile;
    myfile.open(patt + ".dot");
    
    std::string dot = "digraph G {\n graph [splines=true overlap=false]\n node  [shape=ellipse, width=0.3, height=0.3]\n";
    
    for (TNodeEDatNet<TInt, TInt>::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) { 
	dot += " " + std::to_string(NI.GetId()) + " [label=\"" + labels.GetDat(labels.GetKey(labels.GetKeyId(NI.GetId()))).GetCStr() + "\"];\n";
    }
    
    for (TNodeEDatNet<TInt, TInt>::TEdgeI EI = Graph->BegEI(); EI < Graph->EndEI(); EI++) { 
	dot += " " + std::to_string(EI.GetSrcNId()) + " -> " + std::to_string(EI.GetDstNId()) + " [label=\" " + std::to_string(EI.GetDat()) + "\"];\n";
    }
    dot += "}";
    
    myfile << dot;
    myfile.close();
    if(system(("dot -Tpng ./" + patt + ".dot -o ./" + patt + ".png").c_str()) != 0) {
	std::cerr << "System call error" << std::endl;
    }
}
