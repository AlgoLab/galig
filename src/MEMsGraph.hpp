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
    int min_w;
    int exsN;
    lemon::ListDigraph graph;
    lemon::ListDigraph::NodeMap<Mem> nodes_map;
    lemon::ListDigraph::ArcMap<int> edges_map;
    lemon::ListDigraph::Node start;
    lemon::ListDigraph::Node end;
    std::vector<lemon::ListDigraph::Node> starting_nodes;
    std::vector<lemon::ListDigraph::Node> ending_nodes;
    std::list<std::list<Mem> > paths;

    void combine_MEMs(const SplicingGraph&,
                      const std::string&,
                      std::list<Mem>,
                      const int&,
                      lemon::ListDigraph&,
                      lemon::ListDigraph::NodeMap<Mem>&,
                      lemon::ListDigraph::ArcMap<int>&);
    void saveImage(const std::string&,
                   const lemon::ListDigraph&,
                   const lemon::ListDigraph::NodeMap<Mem>&,
                   const lemon::ListDigraph::ArcMap<int>&);
public:
    MemsGraph(const SplicingGraph&,
              const std::string&,
              std::list<Mem>&,
              const int&,
              const int&);
    void build();
    std::pair<int, std::list<std::list<Mem> > > visit();
    
};

#endif
