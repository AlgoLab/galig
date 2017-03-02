#ifndef _FASTA_READER_HPP_
#define _FASTA_READER_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>

class FastaReader {
private:
  std::vector<std::string> descriptions;
  std::vector<std::string> sequences;
public:
  FastaReader(const std::string&);
  int getSize();
  std::pair<std::string, std::string> getEntry(const int&);
};
#endif
