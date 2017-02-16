//=================================
// include guard
#ifndef _MEMSGRAPH_HPP_
#define _MEMSGRAPH_HPP_

//=================================
// included dependencies
#include <string>
#include <unordered_map>
#include <list>

#include "MEMsList.hpp"
#include "SplicingGraph.hpp"

#include <lemon/list_graph.h>

class MemsGraph {
private:
  
  int nodes_index = 0;
  int plen = 0;
  float perc = 0;
public:
  MemsGraph(SplicingGraph &g, MemsList& ml, const int& K, const int& perc);
};

#endif
