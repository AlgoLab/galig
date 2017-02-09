//=================================
// include guard
#ifndef _MEMSGRAPH_HPP_
#define _MEMSGRAPH_HPP_

//=================================
// included dependencies
#include <string>
#include <unordered_map>
#include <list>

#include "utils.hpp"
#include "SplicingGraph.hpp"
#include "MEMsList.hpp"

#include "Snap.h"

class MemsGraph {
private:
    int nodes_index = 0;
    int plen = 0;
    int K = 0;
    float perc = 0;
    PNGraph Graph;
    TIntStrH labels;
    TStrIntH weights;
    std::unordered_map<std::string, int> MemToIndex;
    std::unordered_map<int, Mem> IndexToMem;
    std::vector<std::vector<std::vector<int> > > subpaths;
    std::vector<std::vector<int> > paths;

    TStr toTStr(const std::string& s);
    void addNode(Mem mem);
    void addEdge(Mem mem1, Mem Mem2, int w);
    int getNodeId(const std::string& mem);
    bool isNode(Mem m);
    std::vector<std::vector<int> > rec_visit(const TNGraph::TNodeI node);
public:
    MemsGraph(SplicingGraph& g, MemsList& ml, const int& K, const int& perc);
    void saveImage(const std::string& patt);
    void saveOutput(std::ostream& os, std::string p);
    void visit();
    std::list<std::string> getOutput();
};

#endif
