#include "MEMsGraph.hpp"
#include <time.h>
#include <stdio.h>

MemsGraph::MemsGraph(const std::string& read_,
                     const int& L_,
                     const int& eps_,
                     const int& exsN_,
                     const bool& verbose_) : nodesMap(graph), edgesMap(graph)  {
    read = read_;
    m = read.size();
    L = L_;
    eps = eps_;
    K0 = (eps*m*L)/100;
    K1 = (eps*m*L)/100 - L + 1;
    K2 = (eps*m)/100;
    exsN = exsN_;
    verbose = verbose_;

    start = graph.addNode();
    end = graph.addNode();
    nodesMap[start] = Mem(0,0,0);
    nodesMap[end] = Mem(-1,-1,-1);
}

std::pair<bool, int> MemsGraph::checkMEMs(const SplicingGraph& sg,
                                          const Mem& m1,
                                          const Mem& m2,
                                          const bool& twoPass) {
    int id1 = sg.rank(m1.t-1);
    int id2 = sg.rank(m2.t-1);

    std::string exon1_text = sg.getExon(id1);
    std::string exon2_text = sg.getExon(id2);

    int err = K2;
    bool flag = false;
    if(id1 == id2) { //m1 and m2 in the same exon
        if(m2.p+m2.l>m1.p+m1.l && m1.t<m2.t && m2.t<m1.t+m1.l+K1 && m1.t+m1.l<m2.t+m2.l) {
            int gap_P = m2.p-m1.p-m1.l;
            int gap_E = m2.t-m1.t-m1.l;
            if(gap_P>=0 && gap_E>=0) {
                if(abs(gap_P-gap_E)<=K2) {
                    std::string sub_P = read.substr(m1.p+m1.l-1, m2.p-m1.p-m1.l);
                    std::string sub_E = exon1_text.substr(m1.t+m1.l-sg.select(id1)-1-1, m2.t-m1.t-m1.l);
                    err = e_distance(sub_P, sub_E);
                }
            } else if(gap_P<=0 && gap_E<=0) {
                err = abs(gap_P-gap_E);
            } else {
                err = abs(gap_P) + abs(gap_E);
            }
            if(err < K2) {
                flag = true;
            }
        }
    } else { //m1 and m2 in different exons
        if(sg.contain(id1, id2)) {
            if(!twoPass) {
                if(m2.p+m2.l>m1.p+m1.l) {
                    std::string sub_E1 = exon1_text.substr(m1.t+m1.l-sg.select(id1)-1-1, sg.select(id1+1)+1-m1.t-m1.l);
                    std::string sub_E2 = exon2_text.substr(0,m2.t-sg.select(id2)-1-1);
                    std::string sub_E = sub_E1 + sub_E2;
                    int len_P = m2.p-m1.p-m1.l;
                    std::string sub_P;
                    if(len_P <= 0) {
                        err = abs(len_P);
                    } else {
                        sub_P = read.substr(m1.p+m1.l-1,len_P);
                        err = e_distance(sub_P, sub_E);
                    }
                    if(err < K2) {
                        flag = true;
                    }
                }
            } else {
                if(m2.p>m1.p+m1.l && m2.p+m2.l>m1.p+m1.l && m1.t+m1.l-1==sg.select(id1+1) && m2.t-1-1==sg.select(id2)) {
                    err = 0;
                    flag = true;
                }
            }
        }
    }
    return std::make_pair(flag, err);
}

std::pair<bool, int> MemsGraph::validStart(const SplicingGraph& sg, const Mem& Mem) {
    if(Mem.p <= K0) {
        int err = K2;
        if(Mem.p == 1) {
            err = 0;
        } else {
            int id = sg.rank(Mem.t-1);
            std::string exon_text = sg.getExon(id);
            int l = Mem.p-1;
            std::string sub_P = read.substr(0,l);
            std::string sub_E;
            int exon_pref_len = Mem.t - sg.select(id)-1 - 1;
            if(exon_pref_len < l) {
                int shared_pref_len = l-exon_pref_len;
                std::string exon_pref = exon_text.substr(0, exon_pref_len);
                err = l;
                std::list<int> parents = sg.getParents(id);
                for(std::list<int>::iterator it=parents.begin(); it!=parents.end(); ++it) {
                    /**
                     * We look only at the father of the node,
                     * IF he is long enough, we get its suffix;
                     * ELSE we get all its label (without going further checking all its parents)
                     **/
                    int par = *it;
                    std::string par_text = sg.getExon(par);
                    if(sg.select(par+1)-shared_pref_len-sg.select(par)-1>=0) {
                        sub_E = par_text.substr(sg.select(par+1)-shared_pref_len-sg.select(par)-1, shared_pref_len) + exon_pref;
                    } else {
                        sub_E = par_text + exon_pref;
                    }
                    int curr_err = e_distance(sub_P, sub_E);
                    if(curr_err < err) {
                        err = curr_err;
                    }
                }
            } else {
                sub_E = exon_text.substr(Mem.t-l-sg.select(id)-1-1, l);
                err = e_distance(sub_P, sub_E);
            }
        }
        if(err<K2) {
            return std::make_pair(true, err);
        }
    }
    return std::make_pair(false,K2+1);
}

std::pair<bool, int> MemsGraph::validEnd(const SplicingGraph& sg, const Mem& Mem) {
    if(Mem.p+Mem.l>=m-K0) {
        int err = K2;
        if(Mem.p+Mem.l == m+1) {
            err = 0;
        } else {
            int id = sg.rank(Mem.t-1);
            std::string exon_text = sg.getExon(id);
            int l = m-(Mem.p+Mem.l)+1;
            std::string sub_P = read.substr(Mem.p+Mem.l-1, l);
            std::string sub_E;
            int exon_suff_len = sg.select(id+1) - (Mem.t+Mem.l) + 1;
            if(exon_suff_len < l) {
                std::list<int> sons = sg.getSons(id);
                int shared_suff_len = l-exon_suff_len;
                std::string exon_suff;
                if(exon_suff_len == 0) {
                    exon_suff = "";
                } else {
                    exon_suff = exon_text.substr(Mem.t+Mem.l-sg.select(id)-1-1, exon_suff_len);
                }
                err = l;
                for(std::list<int>::iterator it=sons.begin(); it!=sons.end(); ++it) {
                    int son = *it;
                    std::string son_text = sg.getExon(son);
                    sub_E = exon_suff + son_text.substr(0, shared_suff_len);
                    int curr_err = e_distance(sub_P, sub_E);
                    if(curr_err < err) {
                        err = curr_err;
                    }
                }
            } else {
                sub_E = exon_text.substr(Mem.t+Mem.l-sg.select(id)-1-1, l);
                err = e_distance(sub_P, sub_E);
            }
        }
        if(err<K2) {
            return std::make_pair(true, err);
        }
    }
    return std::make_pair(false,K2+1);
}

void MemsGraph::build(const SplicingGraph& sg,
                      std::list<Mem>& MEMs_,
                      const bool& twoPass) {
    std::vector<std::forward_list<Mem> > MEMs (m+1);
    for(const Mem& m : MEMs_) {
        int p = m.p;
        MEMs[p].push_front(m);
    }
    for(std::vector<std::forward_list<Mem> >::iterator it=MEMs.begin(); it!=MEMs.end(); ++it) {
        for(Mem& m1 : *it) {
            int p1 = m1.p;
            Node node1;
            int id1 = sg.rank(m1.t-1);
            std::string exon1_text = sg.getExon(id1);

            /*********
             * Start *
             *********/
            std::pair<bool, int> start_info = validStart(sg, m1);
            bool start_flag = start_info.first;
            int err = start_info.second;
            if(start_flag) {
                if(err <= K2) {
                    //Possible transitive closures
                    if(m1.isNew) {
                        node1 = graph.addNode();
                        m1.setNode(node1);
                        nodesMap[node1] = m1;
                    } else {
                        node1 = m1.node;
                    }
                    Arc arc = graph.addArc(start,node1);
                    edgesMap[arc] = err;
                } else {
                    continue;
                }
            } else {
                if(m1.isNew) {
                    continue;
                } else {
                    node1 = m1.node;
                }
            }
            /*************
             * Extending *
             *************/
            int p2 = p1+1;
            int max_p = p1+m1.l+K1;
            while((!twoPass && p2<max_p && p2<m) || (twoPass && p2<m)) {
                for(Mem& m2 : MEMs[p2]) {
                    Node node2;
                    std::pair<bool, int> linkage_info = checkMEMs(sg, m1, m2, twoPass);
                    bool flag = linkage_info.first;
                    int err = linkage_info.second;
                    if(flag) {
                        if(m2.isNew) {
                            node2 = graph.addNode();
                            m2.setNode(node2);
                            nodesMap[node2] = m2;
                        } else {
                            node2 = m2.node;
                        }
                        Arc arc = graph.addArc(node1,node2);
                        edgesMap[arc] = err;
                    }
                }
                ++p2;
            }
            /*******
             * End *
             *******/
            std::pair<bool, int> end_info = validEnd(sg, m1);
            bool end_flag = end_info.first;
            err = end_info.second;
            if(end_flag && err <= K2) {
                Arc arc = graph.addArc(node1,end);
                edgesMap[arc] = err;
            }
        }
    }
    /**********************************
     * Transitive closure on end node *
     **********************************/
    std::list<InArc> arcs_D;
    for(InArc XZ (graph, end); XZ!=lemon::INVALID; ++XZ) {
        Node X = graph.source(XZ);
        for(OutArc XY (graph, X); XY!=lemon::INVALID; ++XY) {
            Node Y = graph.target(XY);
            if(graph.id(Y)!=graph.id(end)) {
                for(OutArc YZ (graph, Y); YZ!=lemon::INVALID; ++YZ) {
                    Node Z = graph.target(YZ);
                    if(graph.id(Z)==graph.id(end)) {
                        if(edgesMap[XY]+edgesMap[YZ]<=edgesMap[XZ]) {
                            arcs_D.push_back(XZ);
                        }
                    }
                }
            }
        }
    }
    for(const InArc& a : arcs_D) {
        if(graph.valid(a)) {
            graph.erase(a);
        }
    }
    if(verbose) {
        save("./Graphs/Graph.dot");
    }
}

std::pair<std::pair<bool, int>, std::list<std::pair<bool, std::list<Mem> > > > MemsGraph::visit(const SplicingGraph& sg) {
    bool atLeastOneAnnotaded = false;
    std::list<std::pair<bool, std::list<Mem> > > paths;
    int min_w = K2+1;
    int curr_w;
    for(OutArc arc (graph, start); arc!=lemon::INVALID; ++arc) {
        Node node = graph.target(arc);
        int w0 = edgesMap[arc];
        do {
            lemon::Dijkstra<Graph, lemon::ListDigraph::ArcMap<int> >
                ::SetStandardHeap<FibH>
                ::SetHeap<FibH,FibM>
                ::Create dijkstra (graph, edgesMap);
            FibM heap_cross_ref (graph);
            FibH heap (heap_cross_ref);
            dijkstra.heap(heap, heap_cross_ref);
            dijkstra.run(node,end);
            if(dijkstra.reached(end)) {
                curr_w = w0+dijkstra.dist(end);
                if(curr_w <= min_w) {
                    min_w = curr_w;
                    Path p = dijkstra.path(end);
                    bool annotated = true;
                    int i = 0;
                    std::list<Mem> path;
                    for(Path::ArcIt it(p); it != lemon::INVALID; ++it) {
                        Arc e = it;
                        Node source = graph.source(e);
                        Node target = graph.target(e);
                        Mem m1 = nodesMap[source];
                        if(graph.id(target) != graph.id(end)) {
                            Mem m2 = nodesMap[target];
                            int id1 = sg.rank(m1.t-1);
                            int id2 = sg.rank(m2.t-1);
                            if(id1 != id2 && sg.isNew(id1, id2)) {
                                annotated = false;
                            }
                        }
                        path.push_back(m1);
                        if(i==0) {
                            if(graph.valid(it)) {
                                graph.erase(it);
                            }
                        }
                        ++i;
                    }
                    if(annotated && !atLeastOneAnnotaded) {
                        atLeastOneAnnotaded = true;
                    }
                    paths.push_back(std::make_pair(annotated, path));
                }
            } else {
                curr_w = min_w+1;
            }
        } while(curr_w <= min_w);
        graph.erase(arc);
    }
    return std::make_pair(std::make_pair(atLeastOneAnnotaded, min_w), paths);
}

void MemsGraph::save(const std::string& s) {
    std::ofstream myfile;
    myfile.open(s);

    std::string dot = "digraph G {\n graph [splines=true overlap=false]\n node  [shape=ellipse, width=0.3, height=0.3]\n";
    for (NodeIt n (graph); n != lemon::INVALID; ++n) {
        dot += " " + std::to_string(graph.id(n)) + " [label=\"" + nodesMap[n].toStr() + "\"];\n";
    }
    for(ArcIt a (graph); a != lemon::INVALID; ++a) {
        dot += " " + std::to_string(graph.id(graph.source(a))) + " -> " + std::to_string(graph.id(graph.target(a))) + "[label=\"" + std::to_string(edgesMap[a]) + "\"];\n";
    }
    dot += "}";

    myfile << dot;
    myfile.close();
}
