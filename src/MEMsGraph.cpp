#include "MEMsGraph.hpp"
#include <time.h>
#include <stdio.h>

MemsGraph::MemsGraph(const std::string& read_,
                     const int& L_,
                     const int& eps_,
                     const int& exsN_,
                     const bool& verbose_) : AnnNodesMap(AnnGraph),
                                             AnnEdgesMap(AnnGraph),
                                             NovNodesMap(NovGraph),
                                             NovEdgesMap(NovGraph) {
    read = read_;
    m = read.size();
    L = L_;
    eps = eps_;
    K0 = (eps*m*L)/100;
    K1 = (eps*m*L)/100 - L + 1;
    K2 = (eps*m)/100;
    exsN = exsN_;
    verbose = verbose_;

    AnnStart = AnnGraph.addNode();
    AnnEnd = AnnGraph.addNode();
    NovStart = NovGraph.addNode();
    NovEnd = NovGraph.addNode();

    AnnNodesMap[AnnStart] = Mem(0,0,0);
    AnnNodesMap[AnnEnd] = Mem(-1,-1,-1);
    NovNodesMap[NovStart] = Mem(0,0,0);
    NovNodesMap[NovEnd] = Mem(-1,-1,-1);
}

std::pair<bool, int> MemsGraph::checkMEMs(const SplicingGraph& sg,
                                          const Mem& m1,
                                          const Mem& m2) {
    int id1 = sg.rank(m1.t-1);
    int id2 = sg.rank(m2.t-1);

    std::string exon1_text = sg.getExon(id1);
    std::string exon2_text = sg.getExon(id2);

    int err = -1;
    bool type = true;
    if(id1 == id2) { //m1 and m2 in the same exon
        if(m2.p+m2.l>m1.p+m1.l && m1.t<m2.t && m2.t<m1.t+m1.l+K1 && m1.t+m1.l<m2.t+m2.l) {
            int gapP = m2.p-m1.p-m1.l;
            int gap_E = m2.t-m1.t-m1.l;
            if(gapP>=0 && gap_E>=0) {
                if(gapP == 0) {
                    //Possible intron
                    err = 0;
                    type = false;
                } else if(abs(gapP-gap_E) <= K2) {
                    //Possible SNV
                    std::string sub_P = read.substr(m1.p+m1.l-1, m2.p-m1.p-m1.l);
                    std::string sub_E = exon1_text.substr(m1.t+m1.l-sg.select(id1)-1-1, m2.t-m1.t-m1.l);
                    err = editDistance(sub_P, sub_E);
                    type = true;
                }
            } else if(gapP<=0 && gap_E<=0) {
                err = abs(gapP-gap_E);
                type = true;
            } else {
                err = abs(gapP) + abs(gap_E);
                type = true;
            }
        }
    } else { //m1 and m2 in different exons
        if(sg.contain(id1, id2)) {
            if(m2.p+m2.l>m1.p+m1.l) {
                int gapP = m2.p-m1.p-m1.l;
                int gapE1 = sg.select(id1+1)+1-m1.t-m1.l;
                int gapE2 = m2.t-sg.select(id2)-1-1;
                if(gapP <= 0) {
                    err = abs(gapP);
                    if(!sg.isNew(id1, id2) && gapE1 == 0 && gapE2 == 0) {
                        type = true;
                    }
                    else if(err <= K2) {
                        //Possible Competing
                        type = false;
                    } else
                        err = -1;
                } else {
                    if(gapE1 == 0 && gapE2 == 0) {
                        //Possible insertion
                        err = 0;
                        type = false;
                    } else {
                        if(abs(gapP-(gapE1+gapE2)) <= K2) {
                            //Possible SNV
                            std::string subP = read.substr(m1.p+m1.l-1, gapP);
                            std::string subE1 = exon1_text.substr(m1.t+m1.l-sg.select(id1)-1-1, gapE1);
                            std::string subE2 = exon2_text.substr(0, gapE2);
                            std::string subE = subE1 + subE2;
                            err = editDistance(subP, subE);
                            type = true;
                        }
                    }
                }
            }
        }
    }
    if(err > K2) {
        err = -1;
    }
    return std::make_pair(type, err);
}

std::pair<bool, int> MemsGraph::validStart(const SplicingGraph& sg, const Mem& Mem) {
    if(Mem.p <= K0) {
        int err = K2+1;
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
                    int curr_err = editDistance(sub_P, sub_E);
                    if(curr_err < err) {
                        err = curr_err;
                    }
                }
            } else {
                sub_E = exon_text.substr(Mem.t-l-sg.select(id)-1-1, l);
                err = editDistance(sub_P, sub_E);
            }
        }
        if(err <= K2) {
            return std::make_pair(true, err);
        }
    }
    return std::make_pair(false,K2+1);
}

std::pair<bool, int> MemsGraph::validEnd(const SplicingGraph& sg, const Mem& Mem) {
    if(Mem.p+Mem.l>=m-K0) {
        int err = K2+1;
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
                    int curr_err = editDistance(sub_P, sub_E);
                    if(curr_err < err) {
                        err = curr_err;
                    }
                }
            } else {
                sub_E = exon_text.substr(Mem.t+Mem.l-sg.select(id)-1-1, l);
                err = editDistance(sub_P, sub_E);
            }
        }
        if(err <= K2) {
            return std::make_pair(true, err);
        }
    }
    return std::make_pair(false,K2+1);
}

void MemsGraph::build(const SplicingGraph& sg,
                      std::list<Mem>& MEMs_) {
    std::vector<std::forward_list<Mem> > MEMs (m+1);
    for(const Mem& m : MEMs_) {
        int p = m.p;
        MEMs[p].push_front(m);
    }
    for(std::vector<std::forward_list<Mem> >::iterator it=MEMs.begin(); it!=MEMs.end(); ++it) {
        for(Mem& m1 : *it) {
            int p1 = m1.p;
            Node AnnNode1;
            Node NovNode1;
            int id1 = sg.rank(m1.t-1);
            std::string exon1_text = sg.getExon(id1);

            /*********
             * Start *
             *********/
            std::pair<bool, int> startInfo = validStart(sg, m1);
            bool startFlag = startInfo.first;
            int err = startInfo.second;
            if(startFlag) {
                if(m1.isNew) {
                    AnnNode1 = AnnGraph.addNode();
                    NovNode1 = NovGraph.addNode();
                    m1.setAnnNode(AnnNode1);
                    m1.setNovNode(NovNode1);
                    AnnNodesMap[AnnNode1] = m1;
                    NovNodesMap[NovNode1] = m1;
                } else {
                    AnnNode1 = m1.AnnNode;
                    NovNode1 = m1.NovNode;
                }
                Arc arc = AnnGraph.addArc(AnnStart,AnnNode1);
                AnnEdgesMap[arc] = err;
                arc = NovGraph.addArc(NovStart,NovNode1);
                NovEdgesMap[arc] = err;
            } else {
                if(m1.isNew) {
                    continue;
                } else {
                    AnnNode1 = m1.AnnNode;
                    NovNode1 = m1.NovNode;
                }
            }
            /*************
             * Extending *
             *************/
            bool AnnExt = false;
            bool NovExt = false;
            int p2 = p1+1;
            while(p2 <= m-L+1) {
                for(Mem& m2 : MEMs[p2]) {
                    Node AnnNode2;
                    Node NovNode2;
                    std::pair<bool, int> linkageInfo = checkMEMs(sg, m1, m2);
                    bool flag = linkageInfo.first;
                    int err = linkageInfo.second;
                    if(err>=0) {
                        if(m2.isNew) {
                            AnnNode2 = AnnGraph.addNode();
                            NovNode2 = NovGraph.addNode();
                            m2.setAnnNode(AnnNode2);
                            m2.setNovNode(NovNode2);
                            AnnNodesMap[AnnNode2] = m2;
                            NovNodesMap[NovNode2] = m2;
                        } else {
                            AnnNode2 = m2.AnnNode;
                            NovNode2 = m2.NovNode;
                        }
                        if(flag) {
                            Arc arc = AnnGraph.addArc(AnnNode1,AnnNode2);
                            AnnEdgesMap[arc] = err;
                            arc = NovGraph.addArc(NovNode1,NovNode2);
                            NovEdgesMap[arc] = err;
                            AnnExt = true;
                            NovExt =true;
                        } else {
                            Arc arc = NovGraph.addArc(NovNode1,NovNode2);
                            NovEdgesMap[arc] = err;
                            NovExt =true;
                        }
                    }
                    /**
                    if(flag) {
                        if(m2.isNew) {
                            AnnNode2 = AnnGraph.addNode();
                            NovNode2 = NovGraph.addNode();
                            m2.setNode(node2);
                            nodesMap[node2] = m2;
                        } else {
                            node2 = m2.node;
                        }
                        Arc arc = NovGraph.addArc(node1,node2);
                        edgesMap[arc] = err;
                        extended = true;
                    }
                    **/
                }
                ++p2;
            }
            /*******
             * End *
             *******/
            if(!AnnExt && !NovExt) {
                std::pair<bool, int> endInfo = validEnd(sg, m1);
                bool endFlag = endInfo.first;
                err = endInfo.second;
                if(endFlag) {
                    Arc arc = AnnGraph.addArc(AnnNode1,AnnEnd);
                    AnnEdgesMap[arc] = err;
                    arc = NovGraph.addArc(NovNode1,NovEnd);
                    NovEdgesMap[arc] = err;
                }
            }
        }
    }

    if(verbose) {
        save("Graph.dot");
    }
}

std::list<std::pair<int, std::list<Mem> > > MemsGraph::visit(const SplicingGraph& sg) {
    std::list<std::pair<int, std::list<Mem> > > paths;
    std::list<Mem> AnnPath1;
    std::list<Mem> AnnPath2;
    std::list<Mem> NovPath;
    int AnnW1 = 0;
    int AnnW2 = 0;
    int NovW = 0;

    //Visiting Annotated Graph
    lemon::Dijkstra<Graph, lemon::ListDigraph::ArcMap<int> >
        ::SetStandardHeap<FibH>
        ::SetHeap<FibH,FibM>
        ::Create AnnDijkstra (AnnGraph, AnnEdgesMap);
    FibM AnnHCR (AnnGraph);
    FibH AnnHeap (AnnHCR);
    AnnDijkstra.heap(AnnHeap, AnnHCR);

    AnnDijkstra.run(AnnStart,AnnEnd);
    if(AnnDijkstra.reached(AnnEnd)) {
        AnnW1 = AnnDijkstra.dist(AnnEnd);
        if(AnnW1 <= K2) {
            Path p = AnnDijkstra.path(AnnEnd);
            bool first = true;
            for(Path::ArcIt it(p); it != lemon::INVALID; ++it) {
                if(first) {
                    AnnGraph.erase(it);
                    first = false;
                }
                Arc e = it;
                Node target = AnnGraph.target(e);
                Mem m = AnnNodesMap[target];
                if(NovGraph.id(target) != NovGraph.id(AnnEnd)) {
                    AnnPath1.push_back(m);
                }
            }
            paths.push_back(std::make_pair(AnnW1, AnnPath1));
        }
    }

    AnnDijkstra.run(AnnStart,AnnEnd);
    if(AnnDijkstra.reached(AnnEnd)) {
        AnnW2 = AnnDijkstra.dist(AnnEnd);
        if(AnnW2 <= K2) {
            Path p = AnnDijkstra.path(AnnEnd);
            for(Path::ArcIt it(p); it != lemon::INVALID; ++it) {
                Arc e = it;
                Node target = AnnGraph.target(e);
                Mem m = AnnNodesMap[target];
                if(NovGraph.id(target) != NovGraph.id(AnnEnd)) {
                    AnnPath2.push_back(m);
                }
            }
            paths.push_back(std::make_pair(AnnW2, AnnPath2));
        }
    }

    //Visiting Novel Graph
    lemon::Dijkstra<Graph, lemon::ListDigraph::ArcMap<int> >
        ::SetStandardHeap<FibH>
        ::SetHeap<FibH,FibM>
        ::Create NovDijkstra (NovGraph, NovEdgesMap);
    FibM NovHCR (NovGraph);
    FibH NovHeap (NovHCR);
    NovDijkstra.heap(NovHeap, NovHCR);

    NovDijkstra.run(NovStart,NovEnd);
    if(NovDijkstra.reached(NovEnd)) {
        NovW = NovDijkstra.dist(NovEnd);
        if(NovW < AnnW1) {
            if(NovW <= K2) {
                Path p = NovDijkstra.path(NovEnd);
                for(Path::ArcIt it(p); it != lemon::INVALID; ++it) {
                    Arc e = it;
                    Node target = NovGraph.target(e);
                    Mem m = NovNodesMap[target];
                    if(NovGraph.id(target) != NovGraph.id(NovEnd)) {
                        NovPath.push_back(m);
                    }
                }
            }
            if(paths.size() > 1)
                paths.pop_back();
            paths.push_front(std::make_pair(NovW, NovPath));
        }
    }
    return paths;
}

void MemsGraph::save(const std::string& s) {
    std::ofstream myfile;
    myfile.open(s);

    std::string dot = "digraph G {\n graph [splines=true overlap=false]\n node  [shape=ellipse, width=0.3, height=0.3]\n";
    for (NodeIt n (NovGraph); n != lemon::INVALID; ++n) {
        dot += " " + std::to_string(NovGraph.id(n)) + " [label=\"" + NovNodesMap[n].toStr() + "\"];\n";
    }
    for(ArcIt a (NovGraph); a != lemon::INVALID; ++a) {
        dot += " " + std::to_string(NovGraph.id(NovGraph.source(a))) + " -> " + std::to_string(NovGraph.id(NovGraph.target(a))) + "[label=\"" + std::to_string(NovEdgesMap[a]) + "\"];\n";
    }
    dot += "}";

    myfile << dot;
    myfile.close();
}
