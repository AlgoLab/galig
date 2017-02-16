#ifndef _MEMSGRAPH_HPP_
#define _MEMSGRAPH_HPP_

#include <string>
#include <list>

#include "utils.hpp"
#include "SplicingGraph.hpp"

#include <lemon/list_graph.h>

class MemsGraph {
private:
  
  int nodes_index = 0;
  int plen = 0;
  float perc = 0;
public:
  MemsGraph(SplicingGraph&, std::list<Mem>&, const int&);
};

#endif
