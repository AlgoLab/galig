//=================================
// include guard
#ifndef _SPLICING_GRAPH_HPP_
#define _SPLICING_GRAPH_HPP_

//=================================
// included dependencies
#include <iostream>
#include <fstream>
#include <list>
#include <utility>
#include <string>
#include <map>

#include "sdsl/bit_vectors.hpp"
#include "FastaReader.hpp"
#include "utils.hpp"

struct Feature {
    int i;
    int index;
    std::string seqid;
    std::string type;
    int start;
    int end;
    bool strand; //+: true; -:false
    std::string id;

    Feature() {
        i = 1;
    }

    Feature(std::string line) {
        i = 1;
        std::string token;
        std::string delimiter = "\t";
        std::size_t pos;
        while ((pos = line.find(delimiter)) != std::string::npos) {
            token = line.substr(0, pos);
            line.erase(0, pos + delimiter.length());
            add(token);
        }
        add(line);
    }

    void add(std::string elem) {
        switch(i) {
        case 1:
            seqid = elem;
            break;
        case 3:
            type = elem;
            break;
        case 4:
            start = std::stoi(elem);
            break;
        case 5:
            end = std::stoi(elem);
            break;
        case 7:
            if(elem.compare("+") == 0) {
                strand = true;
            } else {
                strand = false;
            }
            break;
        case 9:
            bool flag = true;
            std::string delimiter = ";";
            std::size_t pos;
            std::string token;
            std::string subtoken;
            while((pos = elem.find(delimiter)) != std::string::npos) {
                token = elem.substr(0, pos);
                if(token.substr(0,2).compare("ID") == 0) {
                    id = token.substr(3,token.size());
                    flag = false;
                    break;
                }
                elem.erase(0, pos + delimiter.length());
            }
            if(flag && elem.substr(0,2).compare("ID") == 0) {
                id = elem.substr(3,elem.size());
                flag = false;
            }
            if(flag) {
                id = ".";
            }
        }
        i++;
    }

    std::string toStr() {
        return seqid + " " +
            type + " " +
            std::to_string(start) + " " +
            std::to_string(end) + " " +
            std::to_string(strand) + " " +
            id;
    }
};

class SplicingGraph {
private:
    std::string T;
    std::vector<std::string> Exons;
    int exsN;
    std::vector<std::vector<int> > edges;
    sdsl::rrr_vector<> bitVector;
    sdsl::rrr_vector<>::select_1_type select_BV;
    sdsl::rrr_vector<>::rank_1_type rank_BV;
    void setupBitVector();
    void save(const std::string);
    void load(const std::string);
public:
    SplicingGraph(const std::string&);
    SplicingGraph(const std::string&, const std::string&, const std::string&);
    std::string getText() const;
    std::string getExon(const int&) const;
    int rank(const int&) const;
    int select(const int&) const;
    bool contain(const std::vector<int>&) const;
    void print() const;
    int getExonsNumber() const;
};

#endif
