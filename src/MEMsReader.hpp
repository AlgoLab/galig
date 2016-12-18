//=================================
// include guard
#ifndef MEMsReader
#define MEMsReader

//=================================
// forward declared dependencies
//class MemsList;

//=================================
// included dependencies
//#include <pair>
#include <string>
#include <forward_list>
#include <vector>
#include <iostream>

#include "MEMsList.hpp"

class MemsReader {
private:
    std::forward_list<std::pair<std::string, MemsList> > patterns;
public:
    MemsReader();
    void addPattern(const std::string& pattern_id, const int& pattern_length, std::forward_list<std::vector<int> > MEMs);
    void getPatterns();
};

#endif
