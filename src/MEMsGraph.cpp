#include "MEMsGraph.hpp"

//ListDigraph::Node x = g.addNode();
//ListDigraph::Arc arc = g.addArc(x,z);

MemsGraph::MemsGraph(const SplicingGraph& sg, const std::string& read, std::list<Mem>& MEMs, const int& L_) {
    m = read.size();
    L = L_;
    int exsN = sg.getExonsNumber();

    lemon::ListDigraph graph;
    lemon::ListDigraph::NodeMap<Mem> nodes_map(graph);
    lemon::ListDigraph::ArcMap<int> edges_map(graph);
    start = graph.addNode();
    nodes_map[start] = Mem();
    end = graph.addNode();
    nodes_map[end] = Mem(-1,-1,-1);

    starting_nodes.resize(exsN);
    ending_nodes.resize(exsN);

    std::vector<std::list<Mem> > divided_MEMs (exsN);
    for(const Mem& m : MEMs) {
        divided_MEMs[sg.rank(m.t)-1].push_back(m);
    }
    int ex_id = 0;
    for(std::list<Mem> mems : divided_MEMs) {
        combine_MEMs(sg, read, mems, ++ex_id, graph, nodes_map, edges_map);
    }

    saveImage("graph", graph, nodes_map, edges_map);
}

void MemsGraph::combine_MEMs(const SplicingGraph& sg,
                             const std::string& read,
                             std::list<Mem> mems,
                             const int& ex_id,
                             lemon::ListDigraph& graph,
                             lemon::ListDigraph::NodeMap<Mem>& nodes_map,
                             lemon::ListDigraph::ArcMap<int>& edges_map) {
    int K0 = L+1; //Inizio/Fine pattern
    int K1 = L+2; //Distanza MEMs
    int K2 = 2; //Confidenza sulla differenza dei gap corti
    int K3 = 3; //Errore permesso nei gap corti
    int K4 = 5; //Errore permesso nei gap lunghi

    std::string exon_text = sg.getExon(ex_id);
    
    lemon::ListDigraph::Node curr_start = graph.addNode();
    nodes_map[curr_start] = Mem(ex_id,ex_id,ex_id);
    starting_nodes[ex_id-1] = curr_start;
    
    std::map<int, std::list<Mem> > ordered_MEMs;
    for(const Mem& m : mems) {
        ordered_MEMs[m.p].push_back(m);
    }
    for(std::map<int, std::list<Mem> >::iterator it=ordered_MEMs.begin(); it!=ordered_MEMs.end(); ++it) {
        it->second.sort(compareMEMs);
    }
    std::map<std::string, lemon::ListDigraph::Node> addedNodes;
    for(std::map<int, std::list<Mem> >::iterator it=ordered_MEMs.begin(); it!=ordered_MEMs.end(); ++it) {
        for(Mem m1 : it->second) {
            lemon::ListDigraph::Node node1;
            try {
                node1 = addedNodes.at(m1.toStr());
            } catch(const std::out_of_range& oor) {
                if(m1.p<=K0 || (m1.p>K0 && (m1.t == sg.select(ex_id)+1+1 || m1.t == sg.select(ex_id)+1+3))) {
                    node1 = graph.addNode();
                    nodes_map[node1] = m1;
                    addedNodes[m1.toStr()] = node1;
                    if(m1.p<=K0) {
                        lemon::ListDigraph::Arc arc = graph.addArc(start,node1);
                        edges_map[arc] = 0;//todo
                    } else {
                        lemon::ListDigraph::Arc arc = graph.addArc(curr_start,node1);
                        edges_map[arc] = 0;//?
                    }
                } else {
                    continue;
                }
            }
            int curr_p = m1.p+1;
            while(curr_p <= m1.p + m1.l + K1) {
                try {
                    for(Mem m2 : ordered_MEMs.at(curr_p)) {
                        if(m1.t<m2.t) {
                            if(m2.t<=m1.t+m1.l+K1 && m1.p+m1.l<m2.p+m2.l && m1.t+m1.l<m2.t+m2.l) {
                                int gap_P = m2.p-m1.p-m1.l;
                                int gap_E = m2.t-m1.t-m1.l;
                                int err;
                                if(gap_P>=0 && gap_E>=0) {
                                    if(abs(gap_P-gap_E)<=K2) {
                                        std::string sub_P = read.substr(m1.p+m1.l-1, m2.p-m1.p-m1.l);
                                        std::string sub_E = exon_text.substr(m1.t+m1.l-sg.select(ex_id)-1-1, m2.t-m1.t-m1.l);
                                        err = e_distance(sub_P, sub_E);
                                    }
                                } else if(gap_P<=0 && gap_E<=0) {
                                    err = abs(gap_P-gap_E);
                                }
                                if(err<=K3) {
                                    lemon::ListDigraph::Node node2;
                                    try {
                                        node2 = addedNodes.at(m1.toStr());
                                    } catch(const std::out_of_range& oor) {
                                        node2 = graph.addNode();
                                        nodes_map[node2] = m2;
                                        addedNodes[m2.toStr()] = node2;
                                    }
                                    lemon::ListDigraph::Arc arc = graph.addArc(node1,node2);
                                    edges_map[arc] = err;
                                }
                            }
                        } else {
                            break;
                        }
                    }
                } catch(const std::out_of_range& oor) {}
                ++curr_p;
            }
        }
    }

    //Adding ENDs
    lemon::ListDigraph::Node curr_end = graph.addNode();
    nodes_map[curr_end] = Mem(ex_id,ex_id,ex_id);
    ending_nodes[ex_id-1] = curr_end;
    for(std::map<std::string, lemon::ListDigraph::Node>::iterator it=addedNodes.begin(); it!=addedNodes.end(); ++it) {
        lemon::ListDigraph::Node node = it->second;
        lemon::ListDigraph::InArcIt in(graph, node);
        Mem m1 = nodes_map[node];
        if(in == lemon::INVALID) {
            if(m1.p+m1.l-1>m-K0 || (m1.p+m1.l-1<=m-K0 && (m1.t+m1.l == sg.select(ex_id+1)+1 || m1.t+m1.l == sg.select(ex_id+1)+1-3))) {
                if(m1.p+m1.l-1>m-K0) {
                    lemon::ListDigraph::Arc arc = graph.addArc(node,end);
                    edges_map[arc] = 0;//todo
                } else {
                    lemon::ListDigraph::Arc arc = graph.addArc(node,curr_end);
                    edges_map[arc] = 0;//?
                }
            }
        }
    }

    //Merging E' -> S' (maybe this can be done in the previous for)
    for(lemon::ListDigraph::InArcIt in (graph, curr_end); in!=lemon::INVALID; ++in) {
        lemon::ListDigraph::Node node1 = graph.target(in);
        Mem m1 = nodes_map[node1];
        for(lemon::ListDigraph::OutArcIt out (graph, curr_start); out!=lemon::INVALID; ++out) {
            lemon::ListDigraph::Node node2 = graph.target(out);
            Mem m2 = nodes_map[node2];

            std::string sub_P = read.substr(m1.p+m1.l-1, m2.p-m1.p-m1.l);
            std::string sub_E = exon_text.substr(m1.t+m1.l-sg.select(ex_id)-1-1, m2.t-m1.t-m1.l);
            int err = e_distance(sub_P, sub_E);
            if(err<=K4) {
                lemon::ListDigraph::Arc arc = graph.addArc(node1,node2);
                edges_map[arc] = err;
                graph.erase(in);
                graph.erase(out);
            }
        }
    }
}

void MemsGraph::saveImage(const std::string& s, const lemon::ListDigraph& graph, const lemon::ListDigraph::NodeMap<Mem>& nodes_map, const lemon::ListDigraph::ArcMap<int>& edges_map) {
    std::ofstream myfile;
    myfile.open(s);
    
    std::string dot = "digraph G {\n graph [splines=true overlap=false]\n node  [shape=ellipse, width=0.3, height=0.3]\n";
    for (lemon::ListDigraph::NodeIt n (graph); n != lemon::INVALID; ++n) {
        dot += " " + std::to_string(graph.id(n)) + " [label=\"" + nodes_map[n].toStr() + "\"];\n";
    }
    for (lemon::ListDigraph::ArcIt a (graph); a != lemon::INVALID; ++a) { 
        dot += " " + std::to_string(graph.id(graph.source(a))) + " -> " + std::to_string(graph.id(graph.target(a))) + "[label=\"" + std::to_string(edges_map[a]) + "\"];\n";
    }
    dot += "}";
    
    myfile << dot;
    myfile.close();
}
