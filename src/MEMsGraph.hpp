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
    int m;
    int L;
    int eps;
    int K0;
    int K1;
    int K2;
    std::string read;
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
              std::list<Mem>&,
              const int&,
              const int&);
    void build(const SplicingGraph&,
               std::list<Mem>&);
    std::pair<int, std::list<std::list<Mem> > > visit();
    void save(const std::string&);
    
};

#endif
