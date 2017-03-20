#include "SplicingGraph.hpp"

std::string getExonID(int s, int e) {
    return std::to_string(s) + ":" + std::to_string(e);
}

SplicingGraph::SplicingGraph(const std::string& f) {
    load(f);
    setupBitVector();
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
    exsN = 0;
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
                    ++exsN;
                }
            }
        }
    }
    addedExons.clear();
    edges.resize(exsN+1);
    parents.resize(exsN+1);
    sons.resize(exsN+1);
    for(int i = 0; i <= exsN; i++) {
        edges[i] = std::vector< int>(exsN+1, 0);  
    }
    T = "|";
    Exons.resize(exsN+1);
    int ex_id = 1;
    std::map<std::string, int> ids_to_index;
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
                    ids_to_index[e.id] = ex_id;
                    std::string curr_ex_string = genomic.substr(e.start-1, e.end-e.start+1);
                    T += curr_ex_string + "|";
                    Exons[ex_id] = curr_ex_string;
                    ++ex_id;
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
                curr_i = ids_to_index[e.id];
                if(last_i != -1) {
                    if(edges[last_i][curr_i] == 0 && last_i != curr_i) {
                        edges[last_i][curr_i] = 1;
                        parents[curr_i].push_back(last_i);
                        sons[last_i].push_back(curr_i);
                    }
                }
                last_i = curr_i;
            }
        }
    }
    int i = 1;
    for(std::pair<int,int> p1 : exs_pos) {
        int j = 1;
        for(std::pair<int,int> p2 : exs_pos) {
            if(p1.second <= p2.first) {
                if(edges[i][j] == 0) {
                    edges[i][j] = 2;
                    parents[j].push_back(i);
                    sons[i].push_back(j);
                }
            }
            j++;
        }
        i++;
    }
    gffFile.close();

    setupBitVector();
    save(fa);
}

void SplicingGraph::setupBitVector() {
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
    return Exons[i];
}

int SplicingGraph::rank(const int& i) const {
    return rank_BV(i);
}

int SplicingGraph::select(const int& i) const {
    return select_BV(i);
}

bool SplicingGraph::contain(const int& x, const int& y) const {
    if(edges[x][y] == 1) {
        return true;
    } else {
        return false;
    }
}

void SplicingGraph::print() const {
    std::cout << T << std::endl;
    for(std::vector<int> v : edges) {
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

void SplicingGraph::save(const std::string path) {
    std::cout << "Saving SG to disk..." << std::endl;
    std::ofstream ofile;
    ofile.open(path + ".sg");
    ofile << T << "\n";
    ofile << edges.size() << "\n";
    for(std::vector<int> v : edges) {
        for(int e : v) {
            ofile << e << " ";
        }
        ofile << "\n";
    }
    ofile.close();
}

void SplicingGraph::load(const std::string path) {
    std::cout << "Loading SG from disk..." << std::endl;
    std::ifstream ifile;
    int c = 0;
    int row = 0;
    std::string line;
    ifile.open(path);
    if(ifile.is_open()) {
        while(getline(ifile,line)) {
            switch(c) {
            case 0: {
                T = line;
                break;
            }
            case 1: {
                exsN = stoi(line);
                edges.resize(exsN);
                for(int i = 0; i < exsN; i++) {
                    edges[i] = std::vector< int>(exsN, 0);
                }
                break;
            }
            default: {
                int column = 0;
                std::string token;
                std::string delimiter = " ";
                std::size_t pos;
                while((pos = line.find(delimiter)) != std::string::npos) {
                    token = line.substr(0, pos);
                    line.erase(0, pos + delimiter.length());
                    edges[row][column] = stoi(token);
                    ++column;
                }
                ++row;
            }
            }
            ++c;
        }
    } else {
        std::cout << "SG not found!!" << std::endl;
        exit(EXIT_FAILURE);
    }
    ifile.close();
}

int SplicingGraph::getExonsNumber() const {
    return exsN;
}
