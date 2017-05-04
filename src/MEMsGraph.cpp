#include "MEMsGraph.hpp"

MemsGraph::MemsGraph(const SplicingGraph& sg,
                     const std::string& read_,
                     std::list<Mem>& MEMs,
                     const int& L_,
                     const int& eps_,
                     const bool& verbose_) : nodes_map(graph), edges_map(graph)  {
    read = read_;
    m = read.size();
    L = L_;
    eps = eps_;
    K0 = (eps*m*L)/100;
    K1 = ((eps*m)/100 - 1)*L + 1;
    K2 = (eps*m)/100;
    exsN = sg.getExonsNumber();
    verbose = verbose_;

    start = graph.addNode();
    end = graph.addNode();
    nodes_map[start] = Mem(0,0,0);
    nodes_map[end] = Mem(-1,-1,-1);

    starting_nodes.resize(exsN);
    ending_nodes.resize(exsN);
}

void MemsGraph::build(const SplicingGraph& sg,
                      std::list<Mem>& MEMs) {
    std::vector<std::list<Mem> > divided_MEMs (exsN);
    for(const Mem& m : MEMs) {
        divided_MEMs[sg.rank(m.t-1)-1].push_back(m);
    }
    int ex_id = 0;
    for(std::list<Mem>& mems : divided_MEMs) {
        combine_MEMs_inside_exon(sg, mems, ++ex_id);
    }
    if(verbose) {
        save("Graphs/Graph1.dot");
    }
    combine_MEMs(sg);
    if(verbose) {
        save("Graphs/Graph2.dot");
    }
    link_start_end(sg);
    if(verbose) {
        save("Graphs/Graph3.dot");
    }
}

void MemsGraph::combine_MEMs_inside_exon(const SplicingGraph& sg,
                                         std::list<Mem> mems,
                                         const int& ex_id) {
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
        for(const Mem& m1 : it->second) {
            /**
             *    Fase 1a - Preprocessing MEM m1
             * Se il MEM m1 è già stato aggiunto, #NON# viene collegato al nodo START
             * (in quanto si creerebbe una chiusura transitiva -> è necessaria una scelta
             * intelligente dei parametri). 
             * Se non è stato già aggiunto, allora lo si aggiunge e lo si collega allo
             * START in qualsiasi caso. In questo secondo caso, se non è "vicino" (K0) 
             * all'inizio del pattern, viene etichettato come non valido.
             **/
            std::map<std::string, lemon::ListDigraph::Node>::iterator it = addedNodes.find(m1.toStr());
            lemon::ListDigraph::Node node1;
            bool start_flag = false;
            if(it != addedNodes.end()) {
                node1 = it->second;
            } else {
                node1 = graph.addNode();
                nodes_map[node1] = m1;
                addedNodes[m1.toStr()] = node1;
                start_flag = true;
            }
            if(start_flag) {
                lemon::ListDigraph::Arc arc = graph.addArc(curr_start,node1);
                edges_map[arc] = 0;
                if(not(m1.p<=K0)) {
                    notValidStart.push_back(node1);
                }
            }
            /**
             *     Fase 1b - Estensione MEM m1
             * Dato m1=(t1,p1,l1), si cerca di estenderlo con qualsiasi altro MEM
             * m2=(t2,p2,l2) tale che p2 in [p1+1, p1+l1+K1]. Si noti che m1 e m2
             * devono essere separati "quasi" (K2) allo stesso modo sia su P che su T.
             * Per stabilire l'errore, bisogna considerare se i due MEMs si overlappano
             * o meno. Se questo errore è "basso" (K3), m2 è aggiunto (se non già aggiunto
             * precedentemente) e si collegano i due MEMs.
             **/
            int curr_p = m1.p+m1.l-K2;
            while(curr_p <= m1.p + m1.l + K1 && curr_p < m && m1.t <= sg.select(ex_id+1)-L+1) {
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
                                        std::string sub_P = read.substr(m1.p+m1.l-1, m2.p-m1.p-m1.l);
                                        if(verbose) {
                                            std::cout << "1" << std::endl;
                                        }
                                        std::string sub_E = exon_text.substr(m1.t+m1.l-sg.select(ex_id)-1-1, m2.t-m1.t-m1.l);
                                        err = e_distance(sub_P, sub_E);
                                        flag = true;
                                    }
                                } else if(gap_P<=0 && gap_E<=0) {
                                    err = abs(gap_P-gap_E);
                                    flag = true;
                                } else {
                                    err = abs(gap_P) + abs(gap_E);
                                    flag = true;
                                }
                                if(flag) {
                                    if(err <= K2) {
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
                                        // x -> y; y -> z; x -> z
                                        for(lemon::ListDigraph::InArcIt xy (graph, node1); xy!=lemon::INVALID; ++xy) {
                                            lemon::ListDigraph::Node x = graph.source(xy);
                                            lemon::ListDigraph::Node y = graph.target(xy);
                                            int w_toY_err = edges_map[xy];
                                            for(lemon::ListDigraph::InArcIt wy (graph, y); wy!=lemon::INVALID; ++wy) {
                                                int w_toY_err_ = edges_map[wy];
                                                if(w_toY_err_<w_toY_err) {
                                                    w_toY_err = w_toY_err_;
                                                }
                                            }
                                            for(lemon::ListDigraph::InArcIt xz (graph, node2); xz!=lemon::INVALID; ++xz) {
                                                lemon::ListDigraph::Node x_ = graph.source(xz);
                                                if(graph.id(x) == graph.id(x_)) {
                                                    if(edges_map[xz]>=w_toY_err + err) { //err = w(y,z)
                                                        if(graph.valid(xz)) {
                                                            graph.erase(xz);
                                                            break;
                                                        }
                                                    } else {
                                                        if(graph.valid(arc)) {
                                                            graph.erase(arc);
                                                            break;
                                                        }
                                                    }
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
        }
    }
    /**
     *     Fase 2 - Aggiunta del nodo END
     * Ogni MEM aggiunto al grafo viene collegato al nodo END solo se
     * non ha figli.
     **/
    lemon::ListDigraph::Node curr_end = graph.addNode();
    nodes_map[curr_end] = Mem(0,-ex_id,0);
    ending_nodes[ex_id-1] = curr_end;
    for(std::map<std::string, lemon::ListDigraph::Node>::iterator it=addedNodes.begin(); it!=addedNodes.end(); ++it) {
        lemon::ListDigraph::Node n1 = it->second;
        lemon::ListDigraph::OutArcIt out1 (graph, n1);
        if(out1 == lemon::INVALID) {
            lemon::ListDigraph::Arc arc = graph.addArc(n1,curr_end);
            edges_map[arc] = 0;
        }
    }
    if(verbose) {
        save("Graphs/g" + std::to_string(ex_id) + "_1.dot");
    }
    /**
     *     Fase 3 - Merging start-end
     **/
    std::list<lemon::ListDigraph::InArcIt> in_arcs1;
    for(lemon::ListDigraph::InArcIt in_arc (graph, curr_end); in_arc!=lemon::INVALID; ++in_arc) {
        lemon::ListDigraph::Node e = graph.source(in_arc);
        Mem Emem = nodes_map[e];
        if(Emem.t+Emem.l < sg.select(ex_id+1)-K0-1) {
            Mem Emem = nodes_map[e];
            for(const lemon::ListDigraph::Node s : notValidStart) {
                Mem Smem = nodes_map[s];
                if(Emem.p < Smem.p && Emem.t < Smem.t && Emem.p + Emem.l <= Smem.p && Emem.t + Emem.l <= Smem.t) {
                    std::string sub_P = read.substr(Emem.p+Emem.l-1, Smem.p-Emem.p-Emem.l);
                    if(verbose) {
                        std::cout << "2" << std::endl;
                    }
                    std::string sub_E = exon_text.substr(Emem.t+Emem.l-sg.select(ex_id)-1-1, Smem.t-Emem.t-Emem.l);
                    int err = e_distance(sub_P, sub_E);
                    if(err <= K2) {
                        lemon::ListDigraph::Arc arc = graph.addArc(e,s);
                        edges_map[arc] = err;
                    }
                }
            }
            if(Emem.p+Emem.l<m-K0) {
                in_arcs1.push_back(in_arc);
            }
        }
    }
    for(const lemon::ListDigraph::InArcIt a : in_arcs1) {
        if(graph.valid(a)) {
            graph.erase(a);
        }
    }
    if(verbose) {
        save("Graphs/g" + std::to_string(ex_id) + "_2.dot");
    }

    /**
     *     Fase 4 - Removing paths
     * Vengono eliminati:
     *   I. i MEMs figli di START che sono interni sia al pattern che all'esone
     *   II. i MEMs padri di END che sono interni sia al pattern che all'esone
     *   III. i MEMs rimasti figli di START da cui non si raggiunge l'END
     **/
    /**
    // I
    std::list<lemon::ListDigraph::OutArcIt> out_arcs1;
    for(lemon::ListDigraph::OutArcIt out_arc (graph, curr_start); out_arc!=lemon::INVALID; ++out_arc) {
        lemon::ListDigraph::Node s = graph.target(out_arc);
        Mem Smem = nodes_map[s];
        if(Smem.p > K0 && Smem.t > sg.select(ex_id) + K0 + 1) {
            out_arcs1.push_back(out_arc);
        }
    }
    for(const lemon::ListDigraph::OutArcIt a : out_arcs1) {
        if(graph.valid(a)) {
            graph.erase(a);
        }
    }

    // II
    std::list<lemon::ListDigraph::InArcIt> in_arcs2;
    for(lemon::ListDigraph::InArcIt in_arc (graph, curr_end); in_arc!=lemon::INVALID; ++in_arc) {
        lemon::ListDigraph::Node e = graph.source(in_arc);
        Mem Emem = nodes_map[e];
        if(Emem.p+Emem.l<m-K0 && Emem.t+Emem.l<sg.select(ex_id+1)-K0) {
            in_arcs2.push_back(in_arc);
        }
    }
    for(const lemon::ListDigraph::InArcIt a : in_arcs2) {
        if(graph.valid(a)) {
            graph.erase(a);
        }
    }

    // III
    std::list<lemon::ListDigraph::OutArcIt> out_arcs2;
    for(lemon::ListDigraph::OutArcIt out_arc (graph, curr_start); out_arc!=lemon::INVALID; ++out_arc) {
        lemon::ListDigraph::Node s = graph.target(out_arc);
        lemon::Bfs<lemon::ListDigraph> bfs(graph);
        bfs.init();
        bfs.addSource(s);
        bfs.start();
        if(!bfs.reached(curr_end)) {
            out_arcs2.push_back(out_arc);
        }
    }
    for(const lemon::ListDigraph::OutArcIt a : out_arcs2) {
        if(graph.valid(a)) {
            graph.erase(a);
        }
    }
    **/
}

void MemsGraph::combine_MEMs(const SplicingGraph& sg) {
    /**
     * I sinks del nodo ex_id vengono combinati con i sources
     * di tutti i suoi figli (sons).
     **/
    int ex_id = 0;
    std::list<lemon::ListDigraph::OutArcIt> starting_arcs_D;
    while(ex_id<exsN) {
        lemon::ListDigraph::Node E_Ex = ending_nodes[ex_id];
        ++ex_id;
        std::string curr_exon_text = sg.getExon(ex_id);
        std::list<int> sons = sg.getSons(ex_id);

        std::list<lemon::ListDigraph::InArcIt> ending_arcs_D;
        for(std::list<int>::iterator it=sons.begin(); it!=sons.end(); ++it) {
            int son = *it;
            lemon::ListDigraph::Node s = starting_nodes[son-1];
            std::string son_text = sg.getExon(son);
            for(lemon::ListDigraph::InArcIt in_arc (graph, E_Ex); in_arc!=lemon::INVALID; ++in_arc) {
                lemon::ListDigraph::Node n1 = graph.source(in_arc);
                Mem m1 = nodes_map[n1];
                if(verbose) {
                    std::cout << "3" << std::endl;
                }
                std::string sub_E_1 = curr_exon_text.substr(m1.t+m1.l-sg.select(ex_id)-1-1, sg.select(ex_id)+1-m1.t-m1.l);
                bool flag = false;
                for(lemon::ListDigraph::OutArcIt out_arc (graph, s); out_arc!=lemon::INVALID; ++out_arc) {
                    lemon::ListDigraph::Node n2 = graph.target(out_arc);
                    Mem m2 = nodes_map[n2];
                    if(verbose) {
                        std::cout << m1.toStr() << " -> " << m2.toStr() << std::endl;
                    }
                    if(m1.p<m2.p && m1.p+m1.l<m2.p+m2.l) { //Il controllo sull'arco è già fatto prendendo solo i sons
                        if(verbose) {
                            std::cout << "4" << std::endl;
                        }
                        std::string sub_E_2 = son_text.substr(0,m2.t-sg.select(son)-1-1);
                        std::string sub_E = sub_E_1 + sub_E_2;
                        int len_P = m2.p-m1.p-m1.l;
                        int err = 0;
                        std::string sub_P;
                        if(len_P == 0) {
                            err = 0;
                        } else if(len_P<0) {
                            err = abs(len_P);
                        } else {
                            sub_P = read.substr(m1.p+m1.l-1,len_P);
                            err = e_distance(sub_P, sub_E);
                        }
                        if(err <= K2) {
                            lemon::ListDigraph::Arc arc = graph.addArc(n1,n2);
                            edges_map[arc] = err;
                            flag = true;
                            starting_arcs_D.push_back(out_arc);
                        }
                    }
                }
                if(flag) {
                    ending_arcs_D.push_back(in_arc);
                }
            }
        }

        for(const lemon::ListDigraph::InArcIt a : ending_arcs_D) {
            if(graph.valid(a)) {
                graph.erase(a);
            }
        }
    }

    for(const lemon::ListDigraph::OutArcIt a : starting_arcs_D) {
        if(graph.valid(a)) {
            graph.erase(a);
        }
    }
}

void MemsGraph::link_start_end(const SplicingGraph& sg) {
    /**
     * Si collegano i vari MEMs ai nodi di START ed END globale.
     * -
     * Tutti i MEMs figli di uno start locale vengono testati
     * per vedere se possono essere figli dello start globale.
     * Vengono testati anche i figli di questi per eventuali
     * sovrapposizioni.
     * -
     * Tutti i MEMs padri di un end locale vengono testati
     * per vedere se possono essere figli dell'end globale.
     * Vengono testati anche i padri di questi per eventuali
     * sovrapposizioni.
     **/
    int ex_id = 0;
    std::list<std::pair<int, lemon::ListDigraph::Node> > starting_arcs_A;
    std::list<std::pair<int, lemon::ListDigraph::Node> > ending_arcs_A;
    std::list<lemon::ListDigraph::OutArcIt> starting_arcs_D;
    std::list<lemon::ListDigraph::InArcIt> ending_arcs_D;
    while(ex_id<exsN) {
        lemon::ListDigraph::Node ex_start = starting_nodes[ex_id];
        lemon::ListDigraph::Node ex_end = ending_nodes[ex_id];
        ++ex_id;
        std::string exon_text = sg.getExon(ex_id);
        std::list<int> parents = sg.getParents(ex_id);
        std::list<int> sons = sg.getSons(ex_id);
        for(lemon::ListDigraph::OutArcIt out_s (graph, ex_start); out_s!=lemon::INVALID; ++out_s) {
            lemon::ListDigraph::Node s = graph.target(out_s);
            Mem Smem = nodes_map[s];
            int w=-1;
            if(Smem.p <= K0+1) {
                int l = Smem.p-1;
                std::string sub_P = read.substr(0,l);
                std::string sub_E;

                if(Smem.p == 1) {
                    w = 0;
                } else {
                    int exon_pref_len = Smem.t - sg.select(ex_id)-1 - 1;
                    if(exon_pref_len < l) {
                        int shared_pref_len = l-exon_pref_len;
                        if(verbose) {
                            std::cout << "5" << std::endl;
                        }
                        std::string exon_pref = exon_text.substr(0, exon_pref_len);
                        w = l;
                        for(const int& p : parents) {
                            //Si guarda solo il padre, se è lungo abbastanza (IF), si prende il suffisso,
                            //altrimenti tutto il testo e basta (senza andare indietro ancora)
                            std::string p_text = sg.getExon(p);
                            if(sg.select(p+1)-shared_pref_len-sg.select(p)-1>=0) {
                                if(verbose) {
                                    std::cout << "6" << std::endl;
                                }
                                sub_E = p_text.substr(sg.select(p+1)-shared_pref_len-sg.select(p)-1, shared_pref_len) + exon_pref;
                            }
                            else {
                                sub_E = p_text + exon_pref;
                            }
                            int curr_w = e_distance(sub_P, sub_E);
                            if(curr_w < w) {
                                w = curr_w;
                            }
                        }
                    }
                    else {
                        if(verbose) {
                            std::cout << "7" << std::endl;
                        }
                        sub_E = exon_text.substr(Smem.t-l-sg.select(ex_id)-1-1, l);
                        w = e_distance(sub_P, sub_E);
                    }
                }
                if(w <= K2) {
                    starting_arcs_A.push_back(std::make_pair(w,s));
                    starting_arcs_D.push_back(out_s);
                }
            }
            //Si rifà la stessa cosa per ogni figlio del nodo
            for(lemon::ListDigraph::OutArcIt out_arc (graph, s); out_arc!=lemon::INVALID; ++out_arc) {
                if(w+edges_map[out_arc]>0) {
                    lemon::ListDigraph::Node s_ = graph.target(out_arc);
                    Mem Smem_ = nodes_map[s_];
                    if(Smem_.t != 0) { // il figlio non è un nodo di end locale
                        int ex_id_ = sg.rank(Smem_.t-1);
                        std::string exon_text_ = sg.getExon(ex_id_);
                        int w_ = w+1;
                        if(Smem_.p <= K0+1) {
                            int l = Smem_.p-1;
                            std::string sub_P = read.substr(0,l);
                            std::string sub_E;
                            if(Smem_.p == 1) {
                                w_ = 0;
                            } else {
                                int exon_pref_len = Smem_.t - sg.select(ex_id_)-1 - 1;
                                if(exon_pref_len < l) {
                                    int shared_pref_len = l-exon_pref_len;
                                    if(verbose) {
                                        std::cout << "5a" << std::endl;
                                    }
                                    std::string exon_pref = exon_text_.substr(0, exon_pref_len);
                                    w_= l;
                                    for(const int& p : parents) {
                                        std::string p_text = sg.getExon(p);
                                        if(sg.select(p+1)-shared_pref_len-sg.select(p)-1>=0) {
                                            if(verbose) {
                                                std::cout << "6a" << std::endl;
                                            }
                                            sub_E = p_text.substr(sg.select(p+1)-shared_pref_len-sg.select(p)-1, shared_pref_len) + exon_pref;
                                        }
                                        else {
                                            sub_E = p_text + exon_pref;
                                        }
                                        int curr_w = e_distance(sub_P, sub_E);
                                        if(curr_w < w_) {
                                            w_ = curr_w;
                                        }
                                    }
                                }
                                else {
                                    if(verbose) {
                                        std::cout << "7a" << std::endl;
                                    }
                                    sub_E = exon_text_.substr(Smem_.t-l-sg.select(ex_id_)-1-1, l);
                                    w_ = e_distance(sub_P, sub_E);
                                }
                            }
                            if(w_ <= K2 && w_<w+edges_map[out_arc]) {
                                starting_arcs_A.push_back(std::make_pair(w_,s_));
                                starting_arcs_D.push_back(out_arc);
                            }
                        }
                    }
                }
            }
        }
        for(lemon::ListDigraph::InArcIt in_e (graph, ex_end); in_e!=lemon::INVALID; ++in_e) {
            lemon::ListDigraph::Node e = graph.source(in_e);
            Mem Emem = nodes_map[e];
            int w=-1;
            if(Emem.p+Emem.l >= m-K0+1) {
                int l = m-(Emem.p+Emem.l)+1;
                if(l == 0) {
                    w = 0;
                } else {
                    std::string sub_P = read.substr(Emem.p+Emem.l-1, l);
                    std::string sub_E;
                    int exon_suff_len = sg.select(ex_id+1) - (Emem.t+Emem.l) + 1;
                    if(exon_suff_len < l) {
                        std::list<int> sons = sg.getSons(ex_id);
                        int shared_suff_len = l-exon_suff_len;
                        std::string exon_suff;
                        if(exon_suff_len == 0) {
                            exon_suff = "";
                        } else {
                            if(verbose) {
                                std::cout << "8" << std::endl;
                            }
                            exon_suff = exon_text.substr(Emem.t+Emem.l-sg.select(ex_id)-1-1, exon_suff_len);
                        }
                        w = l;
                        for(const int& s : sons) {
                            std::string s_text = sg.getExon(s);
                            if(verbose) {
                                std::cout << "9" << std::endl;
                            }
                            sub_E = exon_suff + s_text.substr(0, shared_suff_len);
                            int curr_w = e_distance(sub_P, sub_E);
                            if(curr_w < w) {
                                w = curr_w;
                            }
                        }
                    }
                    else {
                        if(verbose) {
                            std::cout << "10" << std::endl;
                        }
                        sub_E = exon_text.substr(Emem.t+Emem.l-sg.select(ex_id)-1-1, l);
                        w = e_distance(sub_P, sub_E);
                    }
                }
                if(w <= K2) {
                    ending_arcs_A.push_back(std::make_pair(w,e));
                    ending_arcs_D.push_back(in_e);
                }
            }
            //Si rifà la stessa cosa per ogni padre del nodo
            for(lemon::ListDigraph::InArcIt in_arc (graph, e); in_arc!=lemon::INVALID; ++in_arc) {
                if(w+edges_map[in_arc]>0) {
                    lemon::ListDigraph::Node e_ = graph.source(in_arc);
                    Mem Emem_ = nodes_map[e_];
                    if(Emem_.t != 0) { // il padre non è un nodo di start locale
                        int ex_id_ = sg.rank(Emem_.t-1);
                        std::string exon_text_ = sg.getExon(ex_id_);
                        if(Emem_.p+Emem_.l >= m-K0+1) {
                            int l = m-(Emem_.p+Emem_.l)+1;
                            int w_;
                            if(l == 0) {
                                w_ = 0;
                            } else {
                                std::string sub_P = read.substr(Emem_.p+Emem_.l-1, l);
                                std::string sub_E;
                                int exon_suff_len = sg.select(ex_id_+1) - (Emem_.t+Emem_.l) + 1;
                                if(exon_suff_len < l) {
                                    //Se il suffisso di esone rimanente è più corto del suffisso del pattern,
                                    //devo concatenarci uno dei figli
                                    std::list<int> sons = sg.getSons(ex_id_);
                                    int shared_suff_len = l-exon_suff_len;
                                    std::string exon_suff;
                                    if(exon_suff_len == 0) {
                                        exon_suff = "";
                                    } else {
                                        exon_suff = exon_text_.substr(Emem_.t+Emem_.l-sg.select(ex_id_)-1-1, exon_suff_len);
                                    }
                                    w_ = l;
                                    for(const int& s : sons) {
                                        std::string s_text = sg.getExon(s);
                                        sub_E = exon_suff + s_text.substr(0, shared_suff_len);
                                        int curr_w = e_distance(sub_P, sub_E);
                                        if(curr_w < w_) {
                                            w_ = curr_w;
                                        }
                                    }
                                } else {
                                    sub_E = exon_text_.substr(Emem_.t+Emem_.l-sg.select(ex_id_)-1-1, l);
                                    w_ = e_distance(sub_P, sub_E);
                                }
                            }
                            if(w_ <= K2 && w_<w+edges_map[in_arc]) {
                                ending_arcs_A.push_back(std::make_pair(w_,e_));
                            }
                        }
                    }
                }
            }
        }
    }
    for(const lemon::ListDigraph::OutArcIt& a : starting_arcs_D) {
        if(graph.valid(a)) {
            graph.erase(a);
        }
    }
    for(const lemon::ListDigraph::InArcIt& a : ending_arcs_D) {
        if(graph.valid(a)) {
            graph.erase(a);
        }
    }
    for(const std::pair<int, lemon::ListDigraph::Node>& p : starting_arcs_A) {
        lemon::ListDigraph::Arc arc = graph.addArc(start,p.second);
        edges_map[arc] = p.first;
    }
    for(const std::pair<int, lemon::ListDigraph::Node>& p : ending_arcs_A) {
        lemon::ListDigraph::Arc arc = graph.addArc(p.second,end);
        edges_map[arc] = p.first;
    }
}

std::pair<int, std::list<std::list<Mem> > > MemsGraph::visit() {
    std::list<std::list<Mem> > paths;
    int min_w = K2+1;
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
                        if(graph.valid(it)) {
                            graph.erase(it);
                        }
                    } else {
                        lemon::ListDigraph::Arc e = it;
                        lemon::ListDigraph::Node source = graph.source(e);
                        path.push_back(nodes_map[source]);
                    }
                    ++i;
                }
                paths.push_back(path);
            }
        } else {
            curr_w = min_w+1;
        }
    } while(curr_w <= min_w);
    return std::make_pair(min_w, paths);
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
