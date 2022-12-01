#ifndef _SPLICING_GRAPH_HPP_
#define _SPLICING_GRAPH_HPP_

#include <iostream>
#include <fstream>
#include <list>
#include <utility>
#include <string>
#include <map>
#include <unordered_map>

#include "sdsl/bit_vectors.hpp"
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
            //GFF
            //std::string delimiter = ";";
            //std::string string_to_search = "ID=";

            //GTF
            std::string delimiter = "; ";
            std::string string_to_search;
            if(type.compare("gene") == 0) {
                string_to_search = "gene_id \"";
            } else if(type.compare("transcript") == 0) {
                string_to_search = "transcript_id \"";
            } else if(type.compare("exon") == 0) {
                string_to_search = "exon_id \"";
            }
            std::size_t pos;
            std::string token;
            std::string subtoken;
            while((pos = elem.find(delimiter)) != std::string::npos) {
                token = elem.substr(0, pos);
                if(token.substr(0,string_to_search.size()).compare(string_to_search) == 0) {
                    //GTF
                    id = token.substr(string_to_search.size(),token.size()-string_to_search.size()-1);
                    //GFF
                    //id = token.substr(string_to_search.size(),token.size()-string_to_search.size());
                    flag = false;
                    break;
                }
                elem.erase(0, pos + delimiter.length());
            }
            if(flag && elem.substr(0,string_to_search.size()).compare(string_to_search) == 0) {
                //GTF
                id = elem.substr(string_to_search.size(),elem.size()-string_to_search.size()-2);
                //GFF
                //id = elem.substr(string_to_search.size(),elem.size()-string_to_search.size()-1);
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
    std::string reference;
    int refLen;
    bool strand;
    std::vector<std::pair<int, int> > ExonsPos;
    std::vector<std::list<int> > parents;
    std::vector<std::list<int> > sons;
    int exsN;

    std::vector<std::vector<int> > edges;
  
    // sdsl::rrr_vector<> bitVector;
    // sdsl::rrr_vector<>::select_1_type selectBV;
    // sdsl::rrr_vector<>::rank_1_type rankBV;
  
    sdsl::sd_vector<> bitVector;
    sdsl::sd_vector<>::select_1_type selectBV;
    sdsl::sd_vector<>::rank_1_type rankBV;
  
    // sdsl::bit_vector bitVector;
    // sdsl::bit_vector::select_1_type selectBV;
    // sdsl::bit_vector::rank_1_type rankBV;
  
    void setupBitVector();
    void save(const std::string);
    //void load(const std::string);
public:
    //SplicingGraph(const std::string&);
    SplicingGraph(const char*, const std::string&);
    std::string getText() const;
    std::list<int> getParents(const int&) const;
    std::list<int> getSons(const int&) const;
    std::string getExon(const int&) const;
    int rank(const int&) const;
    int select(const int&) const;
    bool contain(const int&, const int&) const;
    bool isNew(const int&, const int&) const;
    void print() const;
    int getExonsNumber() const;
};

#endif
