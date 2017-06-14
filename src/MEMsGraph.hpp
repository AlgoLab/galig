#ifndef _MEMSGRAPH_HPP_
#define _MEMSGRAPH_HPP_

#include <string>
#include <list>
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

    void combine_MEMs_inside_exon(const SplicingGraph&,
                                  std::list<Mem>,
                                  const int&);
    void combine_MEMs(const SplicingGraph&);
    void link_start_end(const SplicingGraph&);
public:
    MemsGraph(const SplicingGraph&,
              const std::string&,
              const int&,
              const int&,
              const bool&);
    void build(const SplicingGraph&,
               std::list<Mem>&);
    std::pair<int, std::list<std::list<Mem> > > visit();
    void save(const std::string&);
};

#endif
