//=================================
// include guard
#ifndef _MEMSGRAPH_HPP_
#define _MEMSGRAPH_HPP_

//=================================
// included dependencies
#include "MEMsList.hpp"
#include "ReferenceGraph.hpp"

#include "Snap.h"

class MemsGraph {
private:
    int nodes_index = 0;
    TPt<TNodeEDatNet<TInt, TInt> > Graph;
    TIntStrH labels;
    std::vector<std::vector<std::vector<int> > > subpaths;

    TStr toTStr(const std::string& s);
    void addNode(const int& exon_index, const std::string& label);
    void addEdge(const int& exon_index_1, const int& exon_index_2, const int& w);
    int getNodeId(const std::string& mem);
    std::vector<std::vector<int> > rec_visit(const TNodeEDatNet<TInt, TInt>::TNodeI node);
public:
    MemsGraph(ReferenceGraph &g, MemsList ml, const int& K);
    void save();
    std::vector<std::vector<int> > visit();
};

#endif
