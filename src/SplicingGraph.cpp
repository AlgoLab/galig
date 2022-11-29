#include "SplicingGraph.hpp"

std::string getExonID(int s, int e) {
    return std::to_string(s) + ":" + std::to_string(e);
}

SplicingGraph::SplicingGraph(const char *genomic,
                             const std::string &gtf) {
    refLen = strlen(genomic);

    std::ifstream gtfFile;

    edges.push_back({0});
    edges.push_back({0});

    parents.push_back({});
    sons.push_back({});
    
    T = "|";

    int exID = 1;
    std::map<std::string, int> id2index;
    int curr_i = 1;
    int last_i = -1;
    
    std::string line;
    gtfFile.open(gtf);
    if(gtfFile.is_open()) {
        while(getline(gtfFile,line)) {
            Feature feat (line);
            if(feat.type.compare("gene") == 0) {
                reference = feat.seqid;
                strand = feat.strand;
            }
            else if(feat.type.compare("transcript") == 0 || feat.type.compare("mRNA") == 0) {
                last_i = -1;
            }
            else if(feat.type.compare("exon") == 0) {
                std::string posID = getExonID(feat.start, feat.end);

                if (id2index.find(posID) == id2index.end()) {
                    id2index[posID] = exID;
                    std::string currExString(genomic + feat.start-1, feat.end-feat.start+1);
                    T += currExString + "|";
                    ExonsPos.push_back(std::make_pair(feat.start, feat.end));
                    parents.push_back({});
                    sons.push_back({});
                    ++exID;
                }

                curr_i = id2index[posID];
                if(last_i != -1) {
                    if(last_i != curr_i) {
                        if(last_i >= (int)edges.size()) {
                            int i = edges.size();
                            while(i<=last_i+1) {
                                edges.push_back({0});
                                ++i;
                            }
                        }

                        if(curr_i >= (int)edges[last_i].size()) {
                            int i = edges[last_i].size();
                            while(i<=curr_i+1) {
                                edges[last_i].push_back(0);
                                ++i;
                            }
                        }
                        //edges[last_i] vector<vector<int>> rightValidVariants(1, {i});
                        edges[last_i][curr_i] = 1;
                        if(curr_i >= (int)parents.size()) {
                            int i = 0;
                            while(i<=curr_i+1) {
                                parents.push_back({});
                                ++i;
                            }
                        }
                        parents[curr_i].push_back(last_i);

                        if(last_i >= (int)sons.size()) {
                            int i = 0;
                            while(i<=last_i+1) {
                                sons.push_back({});
                                ++i;
                            }
                        }
                        sons[last_i].push_back(curr_i);
                    }
                }
                last_i = curr_i;
            }
        }
    }

    // Transitive closure on the graph
    int i = 1;
    for(const std::pair<int,int>& p1 : ExonsPos) {
        if(i>=(int)edges.size()) {
            int i_ = edges.size();
            while(i_<=i+1) {
                edges.push_back({0});
                ++i_;
            }
        }
        int j = 1;
        for(const std::pair<int,int>& p2 : ExonsPos) {
            if(p1.second <= p2.first) {
                if(j>=(int)edges[i].size()) {
                    int j_ = edges[i].size();
                    while(j_<=j+1) {
                        edges[i].push_back(0);
                        ++j_;
                    }
                }

                if(edges[i][j] == 0) {
                    edges[i][j] = 2;
                    parents[j].push_back(i);
                    sons[i].push_back(j);
                }
            }
            ++j;
        }
        ++i;
    }
    gtfFile.close();
    setupBitVector();
    save(gtf);
}

void SplicingGraph::setupBitVector() {
    sdsl::bit_vector BV (T.size(), 0);
    unsigned int i = 0;
    while(i<T.size()) {
        if(T[i] == '|') {
            BV[i] = 1;
        }
        i++;
    }
    bitVector = sdsl::rrr_vector<>(BV);
    selectBV = sdsl::rrr_vector<>::select_1_type(&bitVector);
    rankBV = sdsl::rrr_vector<>::rank_1_type(&bitVector);
}

std::string SplicingGraph::getText() const {
    return T;
}

std::list<int> SplicingGraph::getParents(const int& i) const {
    return parents[i];
}

std::list<int> SplicingGraph::getSons(const int& i) const {
    return sons[i];
}

std::string SplicingGraph::getExon(const int& i) const {
    int s = selectBV(i)+1;
    int e = selectBV(i+1);
    std::string exonText = T.substr(s, e-s);
    return exonText;
}

int SplicingGraph::rank(const int& i) const {
    return rankBV(i);
}

int SplicingGraph::select(const int& i) const {
    return selectBV(i);
}

bool SplicingGraph::contain(const int& x, const int& y) const {
    if(edges[x][y] == 1 || edges[x][y] == 2) {
        return true;
    } else {
        return false;
    }
}

bool SplicingGraph::isNew(const int& x, const int& y) const {
    if(edges[x][y] == 2) {
        return true;
    } else {
        return false;
    }
}

void SplicingGraph::print() const {
    std::cout << T << std::endl;
    std::cout << std::endl;
    unsigned int i = 0;
    while(i<bitVector.size()) {
        std::cout << bitVector[i];
        i++;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    for(std::vector<int> v : edges) {
        for(const int& e : v) {
            std::cout << e << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    for(auto e : parents) {
        for(auto p : e) {
            std::cout << p << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    for(auto e : sons) {
        for(auto s : e) {
            std::cout << s << " ";
        }
        std::cout << std::endl;
    }
}

void SplicingGraph::save(const std::string path) {
    std::ofstream ofile;
    ofile.open(path + ".sg");
    ofile << reference << " " << refLen << " ";
    if(strand)
      ofile << "+";
    else
      ofile << "-";
    ofile << "\n";
    ofile << T << "\n";
    ofile << exsN << "\n";
    for(const std::vector<int>& v : edges) {
        for(const int& e : v) {
            ofile << e << " ";
        }
        ofile << "; ";
    }
    ofile << "\n";
    for(const std::pair<int,int>& p : ExonsPos) {
        ofile << p.first << "," << p.second << " ";
    }
    ofile << "\n";
    // for(const std::string& name : ExonsName) {
    //     ofile << name << " ";
    // }
    ofile << "\n";
    ofile.close();
}

int SplicingGraph::getExonsNumber() const {
    return exsN;
}
