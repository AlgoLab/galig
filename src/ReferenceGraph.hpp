//=================================
// include guard
#ifndef _REFERENCEGRAPH_HPP_
#define _REFERENCEGRAPH_HPP_

//=================================
// included dependencies
#include <string>
#include <vector>

#include "../backwardMEM/libs/sdsl-lite/include/sdsl/bit_vectors.hpp"

class ReferenceGraph {
private:
  std::vector< std::vector< int > > edges;
  sdsl::rrr_vector<> bitVector;
  sdsl::rrr_vector<>::select_1_type select_BV;
  sdsl::rrr_vector<>::rank_1_type rank_BV;

  std::vector<int> extractEdge(std::string line);
  std::vector<int> extractExonsLengths(const std::string& fpath);
  void setupEdges(const std::string& fpath1, const std::string& fpath2, int nex);
  int setupBitVector(const std::string& fpath);
public:
  ReferenceGraph(const std::string& exons_file_path, const std::string& real_edges_file_path, const std::string& added_edges_file_path);
  int rank(const int& i);
  int select(const int& i);
  bool contain(std::vector<int> edge);
};

#endif
