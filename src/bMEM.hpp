//=================================
// include guard
#ifndef _BMEM_HPP_
#define _BMEM_HPP_

//=================================
// included dependencies

#include "sdsl/suffix_trees.hpp"
#include "sdsl/suffix_array_algorithm.hpp"
#include "sdsl/csa_wt.hpp"
//#include "../src/testutils.hpp"
#include "sdsl/io.hpp"
#include "sdsl/construct.hpp"

#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <utility>
#include <iomanip>

#include <getopt.h>
#include <cctype>

#include "utils.hpp"

#ifdef BWTK
typedef sdsl::csa_wt<sdsl::csa_wt<>::wavelet_tree_type, BWTK, 10000> tCSA_WT;
#else
typedef sdsl::csa_wt<sdsl::csa_wt<>::wavelet_tree_type, 4, 10000> tCSA_WT;
#endif

#ifdef LCPCOMPRESS
typedef sdsl::cst_sct3<tCSA_WT, sdsl::lcp_dac<> > tCST;
#else
typedef sdsl::cst_sct3<tCSA_WT, sdsl::lcp_dac<> > tCST;
#endif

typedef tCST::size_type size_type;
typedef tCST::node_type node_type;

struct path_item{
  size_type c, lb, rb, p;
  path_item(size_type c=0,
            size_type lb=0,
            size_type rb=0,
            size_type p=0) : c(c), lb(lb), rb(rb), p(p) {
  }
};

typedef std::list<path_item> path_type;

class BackwardMEM {
private:
  tCST cst;
  std::list<Mem> MEMs;
  tCST::size_type backward_search(const tCST::csa_type&,
                                  const unsigned char&,
                                  tCST::size_type&,
                                  tCST::size_type&);
  void addMEM(const std::string&,
              const size_type&,
              const size_type&,
              const size_type&);
public:
  BackwardMEM(const std::string&, const std::string&);
  std::list<Mem> getMEMs(const std::string&, const unsigned int&);
  void printMEMs();
};

#endif
