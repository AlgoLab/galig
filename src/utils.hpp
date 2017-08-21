#ifndef _UTILS_HPP_
#define _UTILS_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include <lemon/list_graph.h>

struct Mem {
    int t;
    int p;
    int l;
    int k;
    bool isNew;
    lemon::ListDigraph::Node AnnNode;
    lemon::ListDigraph::Node NovNode;

    Mem() {
        t = 0;
        p = 0;
        l = 0;
        k = 0;
        isNew = true;
    }

    Mem(int t_, int p_, int l_) {
        t = t_;
        p = p_;
        l = l_;
        k = 0;
        isNew = true;
    }

    void setAnnNode(lemon::ListDigraph::Node node_) {
        AnnNode = node_;
        isNew = false;
    }

    void setNovNode(lemon::ListDigraph::Node node_) {
        NovNode = node_;
        isNew = false;
    }

    std::string toStr() const {
        return "(" + std::to_string(t) + "," +
            std::to_string(p) + "," +
            std::to_string(l) + ")";
    }
};

int editDistance(const std::string&, const std::string&);
std::string reverseAndComplement(const std::string&);
bool compareMEMs(const Mem&, const Mem&);
bool compareMEMsLength(const Mem&, const Mem&);
#endif
