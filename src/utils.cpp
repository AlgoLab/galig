#include "utils.hpp"

std::string reverse_and_complement(const std::string& s) {
  std::string rs = "";
  for(char c : s) {
    switch(c) {
    case 'A':
      rs += "T";
      break;
    case 'C':
      rs += "G";
      break;
    case 'G':
      rs += "C";
      break;
    case 'T':
      rs += "A";
      break;
    default:
      break;
    }
  }
  std::reverse(rs.begin(), rs.end());
  return rs;
}
