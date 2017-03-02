#ifndef _MEMSGRAPH_HPP_
#define _MEMSGRAPH_HPP_

#include <string>
#include <list>

#include "utils.hpp"
#include "SplicingGraph.hpp"

#include <lemon/list_graph.h>

class MemsGraph {
private:
    int m;
    int L;
    lemon::ListDigraph::Node start;
    lemon::ListDigraph::Node end;
    std::vector<lemon::ListDigraph::Node> starting_nodes;
    std::vector<lemon::ListDigraph::Node> ending_nodes;

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
    MemsGraph(const SplicingGraph&, const std::string&, std::list<Mem>&, const int&);
};

#endif
