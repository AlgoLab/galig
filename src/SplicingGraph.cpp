#include "SplicingGraph.hpp"

std::string getExonID(int s, int e) {
    return std::to_string(s) + ":" + std::to_string(e);
}

SplicingGraph::SplicingGraph(const std::string& fa,
                             const std::string& gtf) {
    std::ifstream fastaFile (fa);
    std::string line;
    std::string head;
    std::string genomic;
    if(fastaFile.is_open()) {
        while(getline(fastaFile,line)) {
            if(line.compare("") != 0) {
                if(line[0] == '>') {
                    head = line.substr(1,line.size()-1);
                } else {
                    genomic += line;
                }
            }
        }
    } else {
        std::cerr << "Genomic file not found!" << std::endl;
    }
    refLen = genomic.size();

    std::ifstream gtfFile;
    std::map<std::string, std::list<std::string> > genes;
    std::map<std::string, std::list<std::string> > transcripts;
    std::map<std::string, Feature> exons;
    std::map<std::string, std::string> addedExons;
    std::string currGene;
    std::string currTr;
    exsN = 0;

    gtfFile.open(gtf);
    if(gtfFile.is_open()) {
        while(getline(gtfFile,line)) {
            Feature f (line);
            if(f.type.compare("gene") == 0) {
                reference = f.seqid;
                currGene = f.id;
                strand = f.strand;
                genes.insert(std::make_pair(currGene, std::list<std::string> ()));
                addedExons.clear();
            }
            else if(f.type.compare("transcript") == 0 || f.type.compare("mRNA") == 0) {
                currTr = f.id;
                genes.at(currGene).push_back(currTr);
                transcripts.insert(std::pair<std::string, std::list<std::string> >(currTr, std::list<std::string> ()));
            }
            else if(f.type.compare("exon") == 0) {
                std::string ex = f.id;
                transcripts.at(currTr).push_back(ex);
                try {
                    addedExons.at(ex);
                } catch(const std::out_of_range& oor) {
                    addedExons[ex] = "";
                    exons.insert(std::make_pair(ex, f));
                    ++exsN;
                }
            }
        }
    }

    addedExons.clear();
    //Questo prima exsN è una stima superiore del vero numero di esoni,
    //dato che lo stesso esone può essere messo in due trascritti con ID diverso.
    //Nella passata successiva, mi baso sulla combinazione start_end
    edges.resize(exsN+1);
    parents.resize(exsN+1);
    sons.resize(exsN+1);
    for(int i = 0; i <= exsN; i++) {
        edges[i] = std::vector< int>(exsN+1, 0);  
    }
    T = "|";

    int exID = 1;
    std::map<std::string, int> id2index;
    int curr_i = 1;
    std::unordered_map<std::string, int> realEdges;
    for(std::map<std::string, std::list<std::string> >::iterator it1=genes.begin(); it1!=genes.end(); ++it1) {
        //Forall gene
        for(std::list<std::string>::iterator it2=it1->second.begin(); it2!=it1->second.end(); ++it2) {
            //Forall transcript in gene
            std::list<int> exonsID;
            int last_i = -1;
            for(std::list<std::string>::iterator it3=transcripts[*it2].begin(); it3!=transcripts[*it2].end(); ++it3) {
                //Forall exon in transcript
                std::string exonID = *it3;
                Feature e = exons[exonID];
                std::string posID = getExonID(e.start, e.end);
                try {
                    addedExons.at(posID);
                    exonsID.push_back(id2index[posID]);
                } catch(const std::out_of_range& oor) {
                    addedExons.insert(std::make_pair(posID, ""));
                    id2index[posID] = exID;
                    exonsID.push_back(exID);
                    std::string currExString = genomic.substr(e.start-1, e.end-e.start+1);
                    T += currExString + "|";
                    ExonsPos.push_back(std::make_pair(e.start, e.end));
                    ExonsName.push_back(e.id);
                    ++exID;
                }
                curr_i = id2index[posID];
                if(last_i != -1) {
                    if(last_i != curr_i && edges[last_i][curr_i] == 0) {
                        realEdges[std::to_string(ExonsPos[last_i-1].second) + std::to_string(ExonsPos[curr_i-1].first)] = 1;
                        edges[last_i][curr_i] = 1;
                        parents[curr_i].push_back(last_i);
                        sons[last_i].push_back(curr_i);
                    }
                }
                last_i = curr_i;
            }
            //Transitive closure on the transcript
            // int i = 0;
            // for(const int& id1 : exonsID) {
            //     int j=0;
            //     for(const int& id2 : exonsID) {
            //         if(i<j && id1 != id2) {
            //             if(edges[id1][id2] == 0) {
            //                 edges[id1][id2] = 2;
            //                 parents[id2].push_back(id1);
            //                 sons[id1].push_back(id2);
            //             }
            //         }
            //         ++j;
            //     }
            //     ++i;
            // }
            exonsID.clear();
        }
    }
    // std::cout << std::endl;
    // for(auto p : ExonsPos) {
    //     std::cout << p.first << " " << p.second << std::endl;
    // }
    // std::cout << std::endl;
    // for(auto it=realEdges.begin(); it!=realEdges.end(); ++it) {
    //     std::cout << it->first << std::endl;
    // }
    // std::cout << std::endl;
    // for(std::vector<int> v : edges) {
    //     for(const int& e : v) {
    //         std::cout << e << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl << std::endl;

    //Transitive closure on the full graph
    int i = 1;
    for(const std::pair<int,int>& p1 : ExonsPos) {
        int start1 = p1.first;
        int end1 = p1.second;
        int j = 1;
        for(const std::pair<int,int>& p2 : ExonsPos) {
            int start2 = p2.first;
            int end2 = p2.second;
            if(end1 <= start2) {
                if(edges[i][j] == 0) {
                    //If the intron end/start of e1/e2 is already an edge (different start/end), it is not a new edge
                    bool isAnnIntron = false; //realEdges.find(std::to_string(end1) + std::to_string(start2)) != realEdges.end();
                    //Here we check for alternative end/start of e1/e2 (but equal start/end)
                    if(!isAnnIntron) {
                        std::list<int> FirstExons;
                        std::list<int> SecondExons;
                        for(const std::pair<int,int>& p3 : ExonsPos) {
                            int start3 = p3.first;
                            int end3 = p3.second;
                            if(start1 == start3 or end1 == end3) {
                                FirstExons.push_back(end3);
                            } else if(start2 == start3 or end2 == end3) {
                                SecondExons.push_back(start3);
                            }
                        }
                        for(const int& e : FirstExons) {
                            for(const int& s : SecondExons) {
                                if(realEdges.find(std::to_string(e) + std::to_string(s)) != realEdges.end()) {
                                    isAnnIntron = true;
                                    break;
                                }
                            }
                        }
                    }
                    if(isAnnIntron) {
                        edges[i][j] = 1;
                    } else {
                        edges[i][j] = 2;
                    }
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
    if(edges[x][y] >= 1) {
        return true;
    } else {
        return false;
    }
}

bool SplicingGraph::isNew(const int& x, const int& y) const {
    if(edges[x][y] > 1) {
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
    for(const std::string& name : ExonsName) {
        ofile << name << " ";
    }
    ofile << "\n";
    ofile.close();
}

int SplicingGraph::getExonsNumber() const {
    return exsN;
}
