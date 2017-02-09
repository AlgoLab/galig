#include "SplicingGraph.hpp"

std::string getExonID(int s, int e) {
    return std::to_string(s) + ":" + std::to_string(e);
}

SplicingGraph::SplicingGraph(const std::string& fa, const std::string& gff) {
    std::string genomic = FastaReader(fa).getEntry(0).second;
    std::ifstream gffFile;
    std::map<std::string, std::list<std::string> > genes;
    std::map<std::string, std::list<std::string> > transcripts;
    std::map<std::string, Feature> exons;
    std::map<std::string, std::string> addedExons;
    std::string line;
    std::string curr_gene = "";
    std::string curr_tr = "";
    int exons_number = 0;
    int e_n = 1;
    gffFile.open(gff);
    if(gffFile.is_open()) {
        while(getline(gffFile,line)) {
            Feature f (line);
            if(f.type.compare("gene") == 0) {
                curr_gene = f.id;
                genes.insert(std::pair<std::string, std::list<std::string> >(curr_gene, std::list<std::string> ()));
                addedExons.clear();
            }
            else if(f.type.compare("transcript") == 0 || f.type.compare("mRNA") == 0) {
                curr_tr = f.id;
                genes.at(curr_gene).push_back(curr_tr);
                transcripts.insert(std::pair<std::string, std::list<std::string> >(curr_tr, std::list<std::string> ()));
            }
            else if(f.type.compare("exon") == 0) {
                std::string pos_ID = getExonID(f.start, f.end);
                std::string ex = f.id;
                if(ex == ".") {
                    try {
                        f.id = addedExons.at(pos_ID);
                    } catch(const std::out_of_range& oor) {
                        f.id = "Exon" + std::to_string(e_n);
                        ++e_n;
                    }
                    ex = f.id;
                }
                transcripts.at(curr_tr).push_back(ex);
                try {
                    addedExons.at(pos_ID);
                } catch(const std::out_of_range& oor) {
                    addedExons[pos_ID] = ex;
                    exons.insert(std::pair<std::string, Feature>(ex, f));
                    ++exons_number;
                }
            }
        }
    }
    addedExons.clear();
    edges.resize(exons_number+1);
    for(int i = 0; i <= exons_number; i++) {
        edges[i] = std::vector< int>(exons_number+1, 0);
    }
    T = "|";
    std::list<std::pair<int, int> > exs_pos;
    int curr_i = 0;
    for(std::map<std::string, std::list<std::string> >::iterator it1=genes.begin(); it1!=genes.end(); ++it1) {
        for(std::list<std::string>::iterator it2=it1->second.begin(); it2!=it1->second.end(); ++it2) {
            int last_i = -1;
            for(std::list<std::string>::iterator it3=transcripts[*it2].begin(); it3!=transcripts[*it2].end(); ++it3) {
                Feature e = exons[*it3];
                try {
                    addedExons.at(e.id);
                } catch(const std::out_of_range& oor) {
                    addedExons.insert(std::pair<std::string, std::string>(e.id, ""));
                    T += genomic.substr(e.start-1, e.end-e.start+1) + "|";
                    curr_i++;
                    exs_pos.push_back(std::pair<int, int> (e.start, e.end));
                    /**
                       if(e.strand) {
                       T += genomic.substr(e.start, e.end-e.start) + "|";
                       }
                       else {
                       T += reverse_and_complement(genomic.substr(e.start, e.end-e.start)) + "|";
                       }
                    **/
                }
                if(last_i != -1) {
                    edges[last_i][curr_i] = 1;
                }
                last_i = curr_i;
            }
        }
        int i = 1;
        for(std::pair<int,int> p1 : exs_pos) {
            int j = 1;
            for(std::pair<int,int> p2 : exs_pos) {
                if(p1.second <= p2.first) {
                    if(edges[i][j] == 0) {
                        edges[i][j] = 2;
                    }
                }
                j++;
            }
            i++;
        }
    }
    gffFile.close();

    sdsl::bit_vector BV (T.size(), 0);
    unsigned int i = 0;
    while(i<T.size()) {
        if(T[i] == '|') { //std::to_string(T[i]).compare("|") == 0) {
            BV[i] = 1;
        }
        i++;
    }
    bitVector = sdsl::rrr_vector<>(BV);
    select_BV = sdsl::rrr_vector<>::select_1_type(&bitVector);
    rank_BV = sdsl::rrr_vector<>::rank_1_type(&bitVector);
}

std::string SplicingGraph::getText() {
    return this->T;
}

int SplicingGraph::rank(const int& i) {
    return rank_BV(i);
}

int SplicingGraph::select(const int& i) {
    return select_BV(i);
}

bool SplicingGraph::contain(const std::vector<int>& edge) {
    if(edges[edge[0]][edge[1]] == 1) {
        return true;
    } else {
        return false;
    }
}

void SplicingGraph::print() {
    std::cout << this->T << std::endl;
    for(std::vector<int> v : this->edges) {
        for(int e : v) {
            std::cout << e << " ";
        }
        std::cout << std::endl;
    }
    unsigned int i = 0;
    while(i<bitVector.size()) {
        std::cout << bitVector[i];
        i++;
    }
    std::cout << std::endl;
}
