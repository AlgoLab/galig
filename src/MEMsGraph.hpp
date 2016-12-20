//=================================
// include guard
#ifndef _MEMSGRAPH_HPP_
#define _MEMSGRAPH_HPP_

//=================================
// included dependencies
#include <unordered_map>

#include "MEMsList.hpp"
#include "ReferenceGraph.hpp"

#include "Snap.h"

class MemsGraph {
private:
    int nodes_index = 0;
    TPt<TNodeEDatNet<TInt, TInt> > Graph;
    std::unordered_map<std::string, int> MEMsToIndex;
    std::vector<std::vector<std::vector<int> > > subpaths;
    std::vector<std::vector<int> > paths;

    TStr toTStr(const std::string& s);
    void addNode(const int& exon_index, const std::string& label);
    void addEdge(const int& exon_index_1, const int& exon_index_2, const int& w);
    int getNodeId(const std::string& mem);
    std::vector<std::vector<int> > rec_visit(const TNodeEDatNet<TInt, TInt>::TNodeI node);
public:
    MemsGraph(ReferenceGraph &g, MemsList& ml, const int& K);
    //void saveImage(const std::string& patt);
    void saveOutput(std::ostream& os);
    void visit();
};

#endif
