#ifndef _UTILS_HPP_
#define _UTILS_HPP_

#include <iostream>
#include <string>
#include <algorithm>

struct Mem {
  int t;
  int p;
  int l;

  Mem(int t_, int p_, int l_) {
    t = t_;
	p = p_;
	l = l_;
  }

  std::string toStr() {
	return "(" + std::to_string(t) + "," +
      std::to_string(p) + "," +
      std::to_string(l) + ")";
  }
};

std::string reverse_and_complement(const std::string&);
#endif
