#ifndef _UTILS_HPP_
#define _UTILS_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

struct Mem {
    int t;
    int p;
    int l;

    Mem() {
        t = 0;
        p = 0;
        l = 0;
    }

    Mem(int t_, int p_, int l_) {
        t = t_;
        p = p_;
        l = l_;
    }

    std::string toStr() const {
        return "(" + std::to_string(t) + "," +
            std::to_string(p) + "," +
            std::to_string(l) + ")";
    }
};

int e_distance(const std::string&, const std::string&);
std::string reverse_and_complement(const std::string&);
bool compareMEMs(const Mem&, const Mem&);
#endif
