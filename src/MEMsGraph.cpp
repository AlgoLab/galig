#include "MEMsGraph.hpp"
#include <time.h>
#include <stdio.h>

MemsGraph::MemsGraph(const std::string& read_,
                     const int& L_,
                     const int& eps_,
                     const int& exsN_,
                     const bool& verbose_) : nodes_map(graph), edges_map(graph)  {
    read = read_;
    m = read.size();
    L = L_;
    eps = eps_;
    K0 = (eps*m*L)/100;
    K1 = ((eps*m)/100 - 1)*L + 1;
    K2 = (eps*m)/100;
    exsN = exsN_;
    verbose = verbose_;

    start = graph.addNode();
    end = graph.addNode();
    nodes_map[start] = Mem(0,0,0);
    nodes_map[end] = Mem(-1,-1,-1);

    starting_nodes.resize(exsN);
    ending_nodes.resize(exsN);
}

std::pair<bool, int> MemsGraph::checkMEMs(const SplicingGraph& sg, const Mem& m1, const Mem& m2) {
    int id1 = sg.rank(m1.t-1);
    int id2 = sg.rank(m2.t-1);

    if(verbose) {
        std::cout << "Exon 1" << std::endl;
    }
    std::string exon1_text = sg.getExon(id1);
    if(verbose) {
        std::cout << "Exon 2" << std::endl;
    }
    std::string exon2_text = sg.getExon(id2);

    int err = K2;
    bool flag = false;
    if(id1 == id2) {
        //m1 and m2 in the same exon
        if(m1.t<m2.t && m2.t<m1.t+m1.l+K1 && m1.t+m1.l<m2.t+m2.l) {
            int gap_P = m2.p-m1.p-m1.l;
            int gap_E = m2.t-m1.t-m1.l;
            if(gap_P>=0 && gap_E>=0) {
                if(abs(gap_P-gap_E)<=K2) {
                    std::string sub_P = read.substr(m1.p+m1.l-1, m2.p-m1.p-m1.l);
                    if(verbose) {
                        std::cout << "1" << std::endl;
                    }
                    std::string sub_E = exon1_text.substr(m1.t+m1.l-sg.select(id1)-1-1, m2.t-m1.t-m1.l);
                    err = e_distance(sub_P, sub_E);
                }
            } else if(gap_P<=0 && gap_E<=0) {
                err = 0; //abs(gap_P-gap_E);
            } else {
                err = abs(gap_P) + abs(gap_E);
            }
            if(err < K2) {
                flag = true;
            }
        }
    } else {
        //m1 and m2 in different exons
        if(sg.contain(id1, id2)) {
            if(m2.p+m2.l>m1.p+m1.l) {
                if(verbose) {
                    std::cout << "2" << std::endl;
                }
                std::string sub_E1 = exon1_text.substr(m1.t+m1.l-sg.select(id1)-1-1, sg.select(id1)+1-m1.t-m1.l);
                std::string sub_E2 = exon2_text.substr(0,m2.t-sg.select(id2)-1-1);
                std::string sub_E = sub_E1 + sub_E2;
                int len_P = m2.p-m1.p-m1.l;
                std::string sub_P;
                if(len_P == 0) {
                    err = 0;
                } else if(len_P<0) {
                    err = 0; //abs(len_P);
                } else {
                    sub_P = read.substr(m1.p+m1.l-1,len_P);
                    err = e_distance(sub_P, sub_E);
                }
                if(err < K2) {
                    flag = true;
                }
            }
        }
    }
    return std::make_pair(flag, err);
}

std::pair<bool, int> MemsGraph::validStart(const SplicingGraph& sg, const Mem& Mem) {
    //if(Mem.p <= L+L/3) {
    if(Mem.p <= K0) {
        int err = K2;
        if(Mem.p == 1) {
            err = 0;
        } else {
            int id = sg.rank(Mem.t-1);
            if(verbose) {
                std::cout << "Exon 3" << std::endl;
            }
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
                    //Si guarda solo il padre, se è lungo abbastanza (IF), si prende il suffisso,
                    //altrimenti tutto il testo e basta (senza andare indietro ancora)
                    int par = *it;
                    if(verbose) {
                        std::cout << "Exon 4" << std::endl;
                    }
                    std::string par_text = sg.getExon(par);
                    if(sg.select(par+1)-shared_pref_len-sg.select(par)-1>=0) {
                        if(verbose) {
                            std::cout << "3" << std::endl;
                        }
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
                if(verbose) {
                    std::cout << "4" << std::endl;
                }
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
    //if(Mem.p+Mem.l>=m-L-L/3) {
    if(Mem.p+Mem.l>=m-K0) {
        int err = K2;
        if(Mem.p+Mem.l == m+1) {
            err = 0;
        } else {
            int id = sg.rank(Mem.t-1);
            if(verbose) {
                std::cout << "Exon 5" << std::endl;
            }
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
                    if(verbose) {
                        std::cout << "5" << std::endl;
                    }
                    exon_suff = exon_text.substr(Mem.t+Mem.l-sg.select(id)-1-1, exon_suff_len);
                }
                err = l;
                for(std::list<int>::iterator it=sons.begin(); it!=sons.end(); ++it) {
                    int son = *it;
                    if(verbose) {
                        std::cout << "Exon 6" << std::endl;
                    }
                    std::string son_text = sg.getExon(son);
                    sub_E = exon_suff + son_text.substr(0, shared_suff_len);
                    int curr_err = e_distance(sub_P, sub_E);
                    if(curr_err < err) {
                        err = curr_err;
                    }
                }
            } else {
                if(verbose) {
                    std::cout << "6" << std::endl;
                }
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

/**************************************************************************************************
 **************************************************************************************************
 ***** GREEDY
 **************************************************************************************************
 **************************************************************************************************/
std::pair<int, std::list<Mem> > MemsGraph::build_greedy(const SplicingGraph& sg,
                                                        std::list<Mem>& MEMs) {
    if(verbose) {
        for(const Mem& m : MEMs) {
            std::cout << m.toStr() << " ";
        }
        std::cout << std::endl;
    }
    
    std::list<Mem> match;
    int error = 0;
    MEMs.sort(compareMEMsLength);

    bool notFound = true;
    while(notFound && !MEMs.empty()) {
        Mem m0 = MEMs.front();
        MEMs.pop_front();
        match.push_front(m0);

        bool start_flag = false;
        bool end_flag = false;

        if(verbose) {
            std::cout << "- Checking " << m0.toStr() << std::endl;
        }
        std::pair<bool, int> start_info = validStart(sg, m0);
        if(!start_info.first) {
            Mem m2 = m0;
            while(!start_flag) {
                for(std::list<Mem>::iterator it=MEMs.begin(); it!=MEMs.end(); ++it) {
                    Mem m1 = *it;
                    if(verbose) {
                        std::cout << "-1 Checking " << m1.toStr() << " -> " << m0.toStr() << std::endl;
                    }
                    std::pair<bool, int> linkage_info = checkMEMs(sg, m1, m2);
                    if(linkage_info.first) {
                        match.push_front(m1);
                        error += linkage_info.second;
                        m2 = m1;
                        std::pair<bool, int> start_info_1 = validStart(sg, m1);
                        if(start_info_1.first) {
                            error+=start_info_1.second;
                            start_flag = true;
                            break;
                        }
                    }
                }
                if(!start_flag) {
                    break;
                }
            }
        } else {
            error+=start_info.second;
            start_flag = true;
        }
        if(verbose) {
            std::cout << "- Start found: " << start_flag << std::endl;
        }
        if(!start_flag) {
            match.clear();
            continue;
        }
        std::pair<bool, int> end_info = validEnd(sg, m0);
        if(!end_info.first) {
            Mem m1 = m0;
            while(!end_flag) {
                for(std::list<Mem>::iterator it=MEMs.begin(); it!=MEMs.end(); ++it) {
                    Mem m2 = *it;
                    if(verbose) {
                        std::cout << "-2 Checking " << m1.toStr() << " -> " << m2.toStr() << std::endl;
                    }
                    std::pair<bool, int> linkage_info = checkMEMs(sg, m1, m2);
                    if(linkage_info.first) {
                        match.push_back(m2);
                        error += linkage_info.second;
                        m1 = m2;
                        std::pair<bool, int> end_info_1 = validEnd(sg, m2);
                        if(end_info_1.first) {
                            error+=end_info_1.second;
                            end_flag = true;
                            break;
                        }
                    }
                }
                if(!end_flag) {
                    break;
                }
            }
        } else {
            error+=end_info.second;
            end_flag = true;
        }
        if(verbose) {
            std::cout << "- End found: " << end_flag <<  std::endl;
        }
        if(start_flag && end_flag) {
            notFound = false;
        } else {
            match.clear();
        }
    }
    if(notFound) {
        error = 2*K2;
    }
    return std::make_pair(error, match);
}

/**************************************************************************************************
 **************************************************************************************************
 ***** EXHAUSTIVE
 **************************************************************************************************
 **************************************************************************************************/
void MemsGraph::build(const SplicingGraph& sg,
                      std::list<Mem>& MEMs_) {
    std::map<int, std::forward_list<Mem> > MEMs;
    for(const Mem& m : MEMs_) {
        int p = m.p;
        std::map<int, std::forward_list<Mem> >::iterator it = MEMs.find(p);
        if (it == MEMs.end()) {
            std::forward_list<Mem> ML;
            MEMs[p] = ML;
        }
        MEMs[p].push_front(m);
    }
    for (std::map<int, std::forward_list<Mem> >::iterator it=MEMs.begin(); it!=MEMs.end(); ++it) {
        int p1 = it->first;
        for(Mem& m1 : it->second) {
            lemon::ListDigraph::Node node1;
            int id1 = sg.rank(m1.t-1);
            std::string exon1_text = sg.getExon(id1);
            // #################################################
            // START
            // #################################################
            std::pair<bool, int> start_info = validStart(sg, m1);
            bool start_flag = start_info.first;
            int err = start_info.second;
            if(start_flag) {
                if(err <= K2) {
                    //Si potrebbero creare chiusure transitive
                    if(m1.isNew) {
                        node1 = graph.addNode();
                        m1.setNode(node1);
                        nodes_map[node1] = m1;
                    } else {
                        node1 = m1.node;
                    }
                    lemon::ListDigraph::Arc arc = graph.addArc(start,node1);
                    edges_map[arc] = err;
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
            // #################################################
            // EXTENDING
            // #################################################
            int p2 = p1+1; //m1.l-1;
            int max_p = p1+m1.l+K1;
            while(p2<max_p && p2<m) {
                for(Mem& m2 : MEMs[p2]) {
                    lemon::ListDigraph::Node node2;
                    if(verbose) {
                        std::cout << "Checking " << m1.toStr() << " -> " << m2.toStr();
                    }
                    std::pair<bool, int> linkage_info = checkMEMs(sg, m1, m2);
                    bool flag = linkage_info.first;
                    int err = linkage_info.second;
                    if(flag) {
                        if(m2.isNew) {
                            node2 = graph.addNode();
                            m2.setNode(node2);
                            nodes_map[node2] = m2;
                        } else {
                            node2 = m2.node;
                        }
                        lemon::ListDigraph::Arc arc = graph.addArc(node1,node2);
                        edges_map[arc] = err;
                    }
                }
                ++p2;
            }
            // #################################################
            // ENDING
            // #################################################
            std::pair<bool, int> end_info = validEnd(sg, m1);
            bool end_flag = end_info.first;
            err = end_info.second;
            if(end_flag && err <= K2) {
                lemon::ListDigraph::Arc arc = graph.addArc(node1,end);
                edges_map[arc] = err;
            }
        }
        ++p1;
    }
    // #################################################
    // TRANSITIVE CLOSURE ON END NODE
    // #################################################
    std::list<lemon::ListDigraph::InArcIt> arcs_D;
    for(lemon::ListDigraph::InArcIt XZ (graph, end); XZ!=lemon::INVALID; ++XZ) {
        lemon::ListDigraph::Node X = graph.source(XZ);
        for(lemon::ListDigraph::OutArcIt XY (graph, X); XY!=lemon::INVALID; ++XY) {
            lemon::ListDigraph::Node Y = graph.target(XY);
            if(graph.id(Y)!=graph.id(end)) {
                for(lemon::ListDigraph::OutArcIt YZ (graph, Y); YZ!=lemon::INVALID; ++YZ) {
                    lemon::ListDigraph::Node Z = graph.target(YZ);
                    if(graph.id(Z)==graph.id(end)) {
                        if(edges_map[XY]+edges_map[YZ]<=edges_map[XZ]) {
                            arcs_D.push_back(XZ);
                        }
                    }
                }
            }
        }
    }
    for(const lemon::ListDigraph::InArcIt& a : arcs_D) {
        if(graph.valid(a)) {
            graph.erase(a);
        }
    }
    
    if(verbose) {
        save("./Graphs/Graph.dot");
    }
}

std::pair<bool, std::pair<int, std::list<std::pair<bool, std::list<Mem> > > > > MemsGraph::visit(const SplicingGraph& sg) {
    bool atLeastOneAnnotaded = false;
    std::list<std::pair<bool, std::list<Mem> > > paths;
    int min_w = K2+1;
    int curr_w;
    for(lemon::ListDigraph::OutArcIt arc (graph, start); arc!=lemon::INVALID; ++arc) {
        lemon::ListDigraph::Node node = graph.target(arc);
        int w0 = edges_map[arc];
        do {
            lemon::Dijkstra<lemon::ListDigraph> dijkstra(graph, edges_map);
            dijkstra.run(node,end);
            if(dijkstra.reached(end)) {
                curr_w = w0+dijkstra.dist(end);
                if(curr_w <= min_w) {
                    min_w = curr_w;
                    lemon::Path<lemon::ListDigraph> p = dijkstra.path(end);
                    bool annotated = true;
                    int i = 0;
                    std::list<Mem> path;
                    for (lemon::Path<lemon::ListDigraph>::ArcIt it(p); it != lemon::INVALID; ++it) {
                        lemon::ListDigraph::Arc e = it;
                        lemon::ListDigraph::Node source = graph.source(e);
                        lemon::ListDigraph::Node target = graph.target(e);
                        Mem m1 = nodes_map[source];
                        if(graph.id(target) != graph.id(end)) {
                            Mem m2 = nodes_map[target];
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

    return std::make_pair(atLeastOneAnnotaded, std::make_pair(min_w, paths));
}


void MemsGraph::save(const std::string& s) {
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
