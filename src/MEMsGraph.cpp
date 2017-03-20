#include "MEMsGraph.hpp"

//ListDigraph::Node x = g.addNode();
//ListDigraph::Arc arc = g.addArc(x,z);
//lemon::ListDigraph::OutArcIt out(graph, node);
//lemon::ListDigraph::InArcIt in(graph, node);

MemsGraph::MemsGraph(const SplicingGraph& sg,
                     const std::string& read,
                     std::list<Mem>& MEMs,
                     const int& L_,
                     const int& eps_) : nodes_map(graph), edges_map(graph)  {
    m = read.size();
    L = L_;
    eps = eps_;
    K0 = (eps*m*L)/100; //Inizio/Fine pattern
    K1 = ((eps*m)/100 - 1)*L + 1; //Distanza MEMs
    K2 = (eps*m)/100; //(m*eps)/100; //Errore permesso nei gap

    exsN = sg.getExonsNumber();

    //lemon::ListDigraph graph;
    //lemon::ListDigraph::NodeMap<Mem> nodes_map (graph);
    //lemon::ListDigraph::ArcMap<int> edges_map (graph);
    start = graph.addNode();
    end = graph.addNode();
    nodes_map[start] = Mem(0,0,0);
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

    saveImage("graphs/graph_pre.dot", graph, nodes_map, edges_map);
    
    ex_id = 0;
    while(ex_id<exsN) {
        std::string curr_exon_text = sg.getExon(ex_id+1);
        std::list<int> sons = sg.getSons(ex_id+1);
        lemon::ListDigraph::Node e = ending_nodes[ex_id];
        std::list<lemon::ListDigraph::InArcIt> in_arcs;
        for(std::list<int>::iterator it=sons.begin(); it!=sons.end(); ++it) {
            int son = *it;
            lemon::ListDigraph::Node s = starting_nodes[son-1];
            std::string son_text = sg.getExon(son);
            for(lemon::ListDigraph::InArcIt in_arc (graph, e); in_arc!=lemon::INVALID; ++in_arc) {
                lemon::ListDigraph::Node n1 = graph.source(in_arc);
                Mem m1 = nodes_map[n1];
                std::string sub_E_1 = curr_exon_text.substr(m1.t+m1.l-sg.select(ex_id+1)-1-1, sg.select(ex_id+1)+1-m1.t-m1.l);
                bool flag = false;
                for(lemon::ListDigraph::OutArcIt out_arc (graph, s); out_arc!=lemon::INVALID; ++out_arc) {
                    lemon::ListDigraph::Node n2 = graph.target(out_arc);
                    Mem m2 = nodes_map[n2];
                    std::string sub_E_2 = son_text.substr(0,m2.t-sg.select(son)-1-1);
                    std::string sub_E = sub_E_1 + sub_E_2;
                    std::string sub_P = read.substr(m1.p+m1.l-1,m2.p-m1.p-m1.l);
                    int err = e_distance(sub_P, sub_E);
                    if(err<K2) {
                        lemon::ListDigraph::Arc arc = graph.addArc(n1,n2);
                        edges_map[arc] = err;
                        flag = true;
                    }
                }
                if(flag) {
                    in_arcs.push_back(in_arc);
                }
            }
        }
        for(const lemon::ListDigraph::InArcIt a : in_arcs) {
            graph.erase(a);
        }
        ++ex_id;
    }
    saveImage("graphs/graph.dot", graph, nodes_map, edges_map);

    min_w = K2;
    std::cout << "Err: " << K2 << std::endl;
    int curr_w;
    do {
        lemon::Dijkstra<lemon::ListDigraph> dijkstra(graph, edges_map);
        dijkstra.run(start,end);
        if(dijkstra.reached(end)) {
            curr_w = dijkstra.dist(end);
            if(curr_w <= min_w) {
                min_w = curr_w;
                lemon::Path<lemon::ListDigraph> p = dijkstra.path(end);
                int i = 0;
                std::list<Mem> path;
                for (lemon::Path<lemon::ListDigraph>::ArcIt it(p); it != lemon::INVALID; ++it) {
                    if(i==0) {
                        graph.erase(it);
                    } else {
                        lemon::ListDigraph::Arc e = it;
                        lemon::ListDigraph::Node source = graph.source(e);
                        path.push_back(nodes_map[source]);
                        //std::cout << nodes_map[source].toStr() << " ";
                    }
                    ++i;
                }
                //std::cout << "\n-" << std::endl;
                paths.push_back(path);
            }
        } else {
            curr_w = min_w+1;
        }
    } while(curr_w <= min_w);
}

std::pair<int, std::list<std::list<Mem> > > MemsGraph::visit() {
    return std::make_pair(min_w, paths);
}

void MemsGraph::combine_MEMs(const SplicingGraph& sg,
                             const std::string& read,
                             std::list<Mem> mems,
                             const int& ex_id,
                             lemon::ListDigraph& graph,
                             lemon::ListDigraph::NodeMap<Mem>& nodes_map,
                             lemon::ListDigraph::ArcMap<int>& edges_map) {
    //std::cout << "##### Exon " << ex_id << std::endl;
    std::string exon_text = sg.getExon(ex_id);
    lemon::ListDigraph::Node curr_start = graph.addNode();
    nodes_map[curr_start] = Mem(0,ex_id,0);
    starting_nodes[ex_id-1] = curr_start;

    std::list<lemon::ListDigraph::Node> notValidStart;
    
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
            /**
             *    Fase 1 - Preprocessing MEM m1
             * Se il MEM m1 è già stato aggiunto, #NON# viene collegato al nodo START. 
             * Se non è stato già aggiunto, allora lo si aggiunge e lo si collega allo
             * START in qualsiasi caso. In questo secondo caso, se non è "vicino" (K0) 
             * all'inizio del pattern, viene etichettato come non valido.
             **/
            std::map<std::string, lemon::ListDigraph::Node>::iterator it = addedNodes.find(m1.toStr());
            lemon::ListDigraph::Node node1;
            bool start_flag = false;
            if(it != addedNodes.end()) {
                node1 = it->second;
                // --- Non viene aggiunto allo start in quanto si creerebbe una "qualche" chiusura transitiva
                //(QUESTO IMPONE UNA SCELTA INTELLIGENTE DEI PARAMETRI)
                //if(m1.p<=K0) { //|| (m1.p>K0 && (m1.t == sg.select(ex_id)+1+1 || m1.t == sg.select(ex_id)+1+3))) {
                //    start_flag = true;
                //}
            } else {
                node1 = graph.addNode();
                nodes_map[node1] = m1;
                addedNodes[m1.toStr()] = node1;
                start_flag = true;
            }
            if(start_flag) {
                //std::cout << "Adding " << m1.toStr() << std::endl;
                lemon::ListDigraph::Arc arc = graph.addArc(curr_start,node1);
                edges_map[arc] = 0;
                if(not(m1.p<=K0)) { //|| (m1.p>K0 && (m1.t == sg.select(ex_id)+1+1 || m1.t == sg.select(ex_id)+1+3)))) {
                    notValidStart.push_back(node1);
                }
            }
            //*******************************************************************************************************
            /**
             *     Fase 2 - Estensione MEM m1
             * Dato m1=(t1,p1,l1), si cerca di estenderlo con qualsiasi altro MEM
             * m2=(t2,p2,l2) tale che p2 in [p1+1, p1+l1+K1]. Si noti che m1 e m2
             * devono essere separati "quasi" (K2) allo stesso modo sia su P che su T.
             * Per stabilire l'errore, bisogna considerare se i due MEMs si overlappano
             * o meno. Se questo errore è "basso" (K3), m2 è aggiunto (se non già aggiunto
             * precedentemente) e si collegano i due MEMs.
             **/
            int curr_p = m1.p+1;
            while(curr_p <= m1.p + m1.l + K1 && curr_p < m) {
                std::map<int, std::list<Mem> >::iterator it = ordered_MEMs.find(curr_p);
                ++curr_p;
                if(it != ordered_MEMs.end()) {
                    std::list<Mem> m2_list = it->second;
                    for(Mem m2 : m2_list) {
                        if(m1.t<m2.t) {
                            if(m2.t<=m1.t+m1.l+K1 && m1.p+m1.l<m2.p+m2.l && m1.t+m1.l<m2.t+m2.l) {
                                int gap_P = m2.p-m1.p-m1.l;
                                int gap_E = m2.t-m1.t-m1.l;
                                int err;
                                bool flag = false;
                                if(gap_P>=0 && gap_E>=0) {
                                    if(abs(gap_P-gap_E)<=K2) {
                                        //std::cout << "1a" << std::endl;
                                        std::string sub_P = read.substr(m1.p+m1.l-1, m2.p-m1.p-m1.l);
                                        //std::cout << "1b" << std::endl;
                                        std::string sub_E = exon_text.substr(m1.t+m1.l-sg.select(ex_id)-1-1, m2.t-m1.t-m1.l);
                                        err = e_distance(sub_P, sub_E);
                                        flag = true;
                                    }
                                } else if(gap_P<=0 && gap_E<=0) {
                                    err = abs(gap_P-gap_E);
                                    flag = true;
                                }
                                if(flag) {
                                    if(err<=K2) {
                                        lemon::ListDigraph::Node node2;
                                        try {
                                            node2 = addedNodes.at(m2.toStr());
                                        } catch(const std::out_of_range& oor) {
                                            node2 = graph.addNode();
                                            nodes_map[node2] = m2;
                                            addedNodes[m2.toStr()] = node2;
                                        }
                                        lemon::ListDigraph::Arc arc = graph.addArc(node1,node2);
                                        edges_map[arc] = err;
                                        // Transitive reduction
                                        for(lemon::ListDigraph::InArcIt in_n1 (graph, node1); in_n1!=lemon::INVALID; ++in_n1) {
                                            lemon::ListDigraph::Node p = graph.source(in_n1);
                                            for(lemon::ListDigraph::InArcIt in_n2 (graph, node2); in_n2!=lemon::INVALID; ++in_n2) {
                                                lemon::ListDigraph::Node p_ = graph.source(in_n2);
                                                if(graph.id(p) == graph.id(p_)) {
                                                    //std::cout << nodes_map[p].toStr() << " -> " << nodes_map[node2].toStr() << std::endl;
                                                    graph.erase(graph.arcFromId(graph.id(in_n2)));
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        } else {
                            break;
                        }
                    }
                }
            }
            //*******************************************************************************************************
        }
    }
    /**
     *     Fase 3 - Aggiunta del nodo END
     * Ogni MEM aggiunto al grafo viene collegato al nodo END solo se
     * non ha figli.
     **/
    std::list<lemon::ListDigraph::Node> nodes;
    //std::cout << "Adding END" << std::endl;
    lemon::ListDigraph::Node curr_end = graph.addNode();
    nodes_map[curr_end] = Mem(0,-ex_id,0);
    ending_nodes[ex_id-1] = curr_end;
    for(std::map<std::string, lemon::ListDigraph::Node>::iterator it=addedNodes.begin(); it!=addedNodes.end(); ++it) {
        lemon::ListDigraph::Node n1 = it->second;
        lemon::ListDigraph::OutArcIt out1 (graph, n1);
        if(out1 == lemon::INVALID) {
            nodes.push_back(n1);
            lemon::ListDigraph::Arc arc = graph.addArc(n1,curr_end);
            edges_map[arc] = 0;
        }
    }

    saveImage("graphs/graph_" + std::to_string(ex_id) + "_1.dot", graph, nodes_map, edges_map);
    
    /**
     *     Fase 4 - Merging start-end
     **/
    //std::cout << "Merging" << std::endl;
    std::list<lemon::ListDigraph::InArcIt> in_arcs;
    for(lemon::ListDigraph::InArcIt in_arc (graph, curr_end); in_arc!=lemon::INVALID; ++in_arc) {
        lemon::ListDigraph::Node e = graph.source(in_arc);
        Mem Emem = nodes_map[e];
        if(Emem.t+Emem.l < sg.select(ex_id+1)-K0-1) {
            Mem Emem = nodes_map[e];
            for(const lemon::ListDigraph::Node s : notValidStart) {
                Mem Smem = nodes_map[s];
                if(Emem.p < Smem.p && Emem.t < Smem.t && Emem.p + Emem.l <= Smem.p && Emem.t + Emem.l <= Smem.t) {
                    //std::cout << "2a" << std::endl;
                    std::string sub_P = read.substr(Emem.p+Emem.l-1, Smem.p-Smem.p-Smem.l);
                    //std::cout << "2b" << std::endl;
                    std::string sub_E = exon_text.substr(Emem.t+Emem.l-sg.select(ex_id)-1-1, Smem.t-Emem.t-Emem.l);
                    int err = e_distance(sub_P, sub_E);
                    if(err<=K2) {
                        lemon::ListDigraph::Arc arc = graph.addArc(e,s);
                        edges_map[arc] = err;
                    }
                }
            }
            if(Emem.p+Emem.l<m-K0) {
                in_arcs.push_back(in_arc);
            }
        }
    }
    for(const lemon::ListDigraph::InArcIt a : in_arcs) {
        graph.erase(a);
    }
    in_arcs.clear();

    saveImage("graphs/graph_" + std::to_string(ex_id) + "_2.dot", graph, nodes_map, edges_map);

    /**
     *     Fase 5 - Removing paths
     **/
    //std::cout << "Removing Source (I)" << std::endl;
    std::list<lemon::ListDigraph::OutArcIt> out_arcs;
    for(lemon::ListDigraph::OutArcIt out_arc (graph, curr_start); out_arc!=lemon::INVALID; ++out_arc) {
        lemon::ListDigraph::Node s = graph.target(out_arc);
        Mem Smem = nodes_map[s];
        //std::cout << Smem.toStr() << std::endl;
        if(Smem.p > K0 && Smem.t > sg.select(ex_id) + K0 + 1) {
            //std::cout << "Erasing source (int)" << Smem.toStr() << std::endl;
            out_arcs.push_back(out_arc);
        }
    }
    for(const lemon::ListDigraph::OutArcIt a : out_arcs) {
        graph.erase(a);
    }
    out_arcs.clear();
    //std::cout << "Removing Sink (I)" << std::endl;
    for(lemon::ListDigraph::InArcIt in_arc (graph, curr_end); in_arc!=lemon::INVALID; ++in_arc) {
        lemon::ListDigraph::Node e = graph.source(in_arc);
        Mem Emem = nodes_map[e];
        //std::cout << Emem.toStr() << std::endl;
        if(Emem.p+Emem.l<m-K0 && Emem.t+Emem.l<sg.select(ex_id+1)-K0) {
            //std::cout << "Erasing sink (int)" << Emem.toStr() << std::endl;
            in_arcs.push_back(in_arc);
        }
    }
    for(const lemon::ListDigraph::InArcIt a : in_arcs) {
        graph.erase(a);
    }
    in_arcs.clear();
    //std::cout << "Removing Source (BFS)" << std::endl;
    for(lemon::ListDigraph::OutArcIt out_arc (graph, curr_start); out_arc!=lemon::INVALID; ++out_arc) {
        lemon::ListDigraph::Node s = graph.target(out_arc);
        //Mem Smem = nodes_map[s];
        //std::cout << Smem.toStr() << std::endl;
        lemon::Bfs<lemon::ListDigraph> bfs(graph);
        bfs.init();
        bfs.addSource(s);
        bfs.start();
        if(!bfs.reached(curr_end)) {
            //std::cout << "Erasing source (bfs)" << Smem.toStr() << std::endl;
            out_arcs.push_back(out_arc);
        }
    }

    for(const lemon::ListDigraph::OutArcIt a : out_arcs) {
        graph.erase(a);
    }
    out_arcs.clear();
    
    saveImage("graphs/graph_" + std::to_string(ex_id) + "_3.dot", graph, nodes_map, edges_map);

    /**
     *     Fase 6
     **/
    //std::cout << "\tMerging start..." << std::endl;
    for(lemon::ListDigraph::OutArcIt out_s (graph, curr_start); out_s!=lemon::INVALID; ++out_s) {
        lemon::ListDigraph::Node s = graph.target(out_s);
        Mem Smem = nodes_map[s];
        //std::cout << Smem.toStr() << std::endl;
        if(Smem.p <= K0+1) {
            lemon::ListDigraph::Arc arc = graph.addArc(start,s);
            int l = Smem.p-1;
            //std::cout << l << std::endl;
            //std::cout << "3a" << std::endl;
            std::string sub_P = read.substr(0,l);
            //std::cout << sub_P << std::endl;
            std::string sub_E;
            int w;
            if(Smem.p == 1) {
                w = 0;
            } else {
                int exon_pref_len = Smem.t - sg.select(ex_id)-1 - 1;
                if(exon_pref_len < l) {
                    //std::cout << "." << std::endl;
                    std::list<int> parents = sg.getParents(ex_id);
                    int shared_pref_len = l-exon_pref_len;
                    //std::cout << shared_pref_len << std::endl;
                    //std::cout << "3b" << std::endl;
                    std::string exon_pref = exon_text.substr(0, exon_pref_len);
                    //std::cout << exon_pref << std::endl;
                    w = l;
                    for(const int& p : parents) {
                        std::string p_text = sg.getExon(p);
                        //std::cout << "3c" << std::endl;
                        if(sg.select(p+1)-shared_pref_len-sg.select(p)-1>=0) {
                            sub_E = p_text.substr(sg.select(p+1)-shared_pref_len-sg.select(p)-1, shared_pref_len) + exon_pref;
                        }
                        else {
                            sub_E = p_text + exon_pref;
                        }
                        //std::cout << "-" << std::endl;
                        int curr_w = e_distance(sub_P, sub_E);
                        if(curr_w < w) {
                            w = curr_w;
                        }
                    }
                }
                else {
                    //std::cout << "3d" << std::endl;
                    sub_E = exon_text.substr(Smem.t-l-sg.select(ex_id)-1, l);
                    w = e_distance(sub_P, sub_E);
                }
                //std::cout << sub_E << std::endl;
            }
            //std::cout << w << std::endl;
            edges_map[arc] = w;
            out_arcs.push_back(out_s);
        }
    }

    //std::cout << "\tMerging ends..." << std::endl;
    for(lemon::ListDigraph::InArcIt in_e (graph, curr_end); in_e!=lemon::INVALID; ++in_e) {
        lemon::ListDigraph::Node e = graph.source(in_e);
        Mem Emem = nodes_map[e];
        //std::cout << Emem.toStr() << std::endl;
        if(Emem.p+Emem.l >= m-K0+1) {
            lemon::ListDigraph::Arc arc = graph.addArc(e,end);
            int l = m-(Emem.p+Emem.l)+1;
            int w;
            if(l == 0) {
                w = 0;
            } else {
                //std::cout << l << std::endl;
                //std::cout << "4a" << std::endl;
                std::string sub_P = read.substr(Emem.p+Emem.l-1, l);
                //std::cout << sub_P << std::endl;
                std::string sub_E;
                int exon_suff_len = sg.select(ex_id+1) - (Emem.t+Emem.l) + 1;
                //std::cout << exon_suff_len << std::endl;
                if(exon_suff_len < l) {
                    std::list<int> sons = sg.getSons(ex_id);
                    int shared_suff_len = l-exon_suff_len;
                    //std::cout << shared_suff_len << std::endl;
                    std::string exon_suff;
                    if(exon_suff_len == 0) {
                        exon_suff = "";
                    } else {
                        //std::cout << "4b" << std::endl;
                        exon_suff = exon_text.substr(Emem.t+Emem.l-sg.select(ex_id)-1-1, exon_suff_len);
                    }
                    //std::cout << "-> " << exon_suff << std::endl;
                    w = l;
                    for(const int& s : sons) {
                        std::string s_text = sg.getExon(s);
                        //std::cout << "4c" << std::endl;
                        sub_E = exon_suff + s_text.substr(0, shared_suff_len);
                        //std::cout << s_text.substr(0, shared_suff_len) << std::endl;
                        int curr_w = e_distance(sub_P, sub_E);
                        if(curr_w < w) {
                            w = curr_w;
                        }
                    }
                }
                else {
                    //std::cout << "4d" << std::endl;
                    sub_E = exon_text.substr(Emem.t+Emem.l-sg.select(ex_id)-1-1, l);
                    w = e_distance(sub_P, sub_E);
                }
                //std::cout << sub_E << std::endl;
            }
            //std::cout << w << std::endl;
            edges_map[arc] = w;
            in_arcs.push_back(in_e);
        }
    }
    for(const lemon::ListDigraph::OutArcIt a : out_arcs) {
        graph.erase(a);
    }
    saveImage("graphs/graph_" + std::to_string(ex_id) + "_4.dot", graph, nodes_map, edges_map);
    for(const lemon::ListDigraph::InArcIt a : in_arcs) {
        graph.erase(a);
    }
    saveImage("graphs/graph_" + std::to_string(ex_id) + "_5.dot", graph, nodes_map, edges_map);
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
