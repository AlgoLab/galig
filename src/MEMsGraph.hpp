#ifndef _MEMSGRAPH_HPP_
#define _MEMSGRAPH_HPP_

#include <string>
#include <list>
#include <forward_list>
#include <utility>

#include "utils.hpp"
#include "SplicingGraph.hpp"

#include <lemon/list_graph.h>
#include <lemon/bfs.h>
#include <lemon/dijkstra.h>

class MemsGraph {
private:
    /**
     * L: L parameter
     * eps: epsilon parameter
     * read: read
     * m: read length
     * K0: start/end proximity
     * K1: gap length
     * K2: errors in gap
     * exsN: exons number
     * graph: MEMsGraph
     * nodes_map: node -> mem
     * edges_map: edge -> weigth
     * start: global starting node
     * end: global ending node
     * starting_nodes: local starting nodes list (exon)
     * ending_nodes: local ending nodes list (exon)
     **/
    int L;
    int eps;
    std::string read;
    bool verbose;
    int m;
    int K0;
    int K1;
    int K2;
    int exsN;
    lemon::ListDigraph graph;
    lemon::ListDigraph::NodeMap<Mem> nodes_map;
    lemon::ListDigraph::ArcMap<int> edges_map;
    lemon::ListDigraph::Node start;
    lemon::ListDigraph::Node end;
    std::vector<lemon::ListDigraph::Node> starting_nodes;
    std::vector<lemon::ListDigraph::Node> ending_nodes;
    std::pair<bool, int> checkMEMs(const SplicingGraph&, const Mem&, const Mem&);
    std::pair<bool, int> validStart(const SplicingGraph&, const Mem&);
    std::pair<bool, int> validEnd(const SplicingGraph&, const Mem&);
public:
    MemsGraph(const std::string&,
              const int&,
              const int&,
              const int&,
              const bool&);
    void build(const SplicingGraph&,
               std::list<Mem>&);
    std::pair<int, std::list<Mem> > build_greedy(const SplicingGraph&,
                                                 std::list<Mem>&);
    std::pair<bool, std::pair<int, std::list<std::pair<bool, std::list<Mem> > > > > visit(const SplicingGraph&);
    void save(const std::string&);
};

#endif
