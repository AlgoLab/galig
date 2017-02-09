//=================================
// include guard
#ifndef _MEMSLIST_HPP_
#define _MEMSLIST_HPP_

//=================================
// forward declared dependencies
//class foo;

//=================================
// included dependencies
#include <iostream>
#include <string>
#include <vector>
#include <forward_list>

#include "utils.hpp"

class MemsList {
private:
    int length;
    std::vector<std::forward_list<Mem> > mems;
public:
    MemsList(const int&);
    void addMem(const Mem&);
    std::forward_list<Mem> getMems(const int&);
    int getLength();
    void print();
};

#endif
