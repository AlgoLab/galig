//=================================
// include guard
#ifndef _MEMSREADER_HPP_
#define _MEMSREADER_HPP_

//=================================
// included dependencies
#include <string>
#include <forward_list>
#include <vector>
#include <iostream>
#include <fstream>

#include "MEMsList.hpp"

class MemsReader {
private:
  std::ifstream memsFile;
  std::forward_list<std::pair<std::string, MemsList> > patterns;
  std::vector<int> extractMEM(std::string line);
  int extractLength(std::string line);
public:
  MemsReader(const std::string& fpath);
  void addPattern(const std::string& pattern_id, const int& pattern_length, std::forward_list<std::vector<int> > MEMs);
  void readMEMsFile();
  void print();
  bool hasPattern();
  std::pair<std::string, MemsList> popPattern();
  int numPatterns();
};

#endif
