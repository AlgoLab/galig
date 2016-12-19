//=================================
// include guard
#ifndef _REFERENCEGRAPH_HPP_
#define _REFERENCEGRAPH_HPP_

//=================================
// included dependencies
#include <string>
#include <vector>

#include "sdsl/bit_vectors.hpp"
#include "Snap.h"

class ReferenceGraph {
private:
    std::vector<std::vector<int> > edges;
    sdsl::rrr_vector<> bitVector;
    
    std::vector<int> extractEdge(std::string line);
    std::vector<int> extractExonsLengths(const std::string& fpath);
    void setupEdges(const std::string& fpath);
    void setupBitVector(const std::string& fpath);
public:
    ReferenceGraph(const std::string& exons_file_path, const std::string& edges_file_path);
};

#endif
