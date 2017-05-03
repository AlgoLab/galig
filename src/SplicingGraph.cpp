#include "SplicingGraph.hpp"

std::string getExonID(int s, int e) {
    return std::to_string(s) + ":" + std::to_string(e);
}

SplicingGraph::SplicingGraph(const std::string& f) {
    load(f);
    setupBitVector();
}

SplicingGraph::SplicingGraph(const std::string& fa,
                             const std::string& gtf) {
    std::string genomic = FastaReader(fa).getEntry(0).second;
    ref_length = genomic.size();

    std::ifstream gtfFile;
    std::map<std::string, std::list<std::string> > genes;
    std::map<std::string, std::list<std::string> > transcripts;
    std::map<std::string, Feature> exons;
    std::map<std::string, std::string> addedExons;
    std::string line;
    std::string curr_gene = "";
    std::string curr_tr = "";
    exsN = 0;
    //int e_n = 1;
    gtfFile.open(gtf);
    if(gtfFile.is_open()) {
        while(getline(gtfFile,line)) {
            Feature f (line);
            if(f.type.compare("gene") == 0) {
                reference = f.seqid;
                curr_gene = f.id;
                genes.insert(std::make_pair(curr_gene, std::list<std::string> ()));
                addedExons.clear();
            }
            else if(f.type.compare("transcript") == 0 || f.type.compare("mRNA") == 0) {
                curr_tr = f.id;
                genes.at(curr_gene).push_back(curr_tr);
                transcripts.insert(std::pair<std::string, std::list<std::string> >(curr_tr, std::list<std::string> ()));
            }
            else if(f.type.compare("exon") == 0) {
                std::string ex = f.id;
                transcripts.at(curr_tr).push_back(ex);
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

    Exons.resize(exsN+1);
    int ex_id = 1;
    std::map<std::string, int> ids_to_index;
    std::list<std::pair<int, int> > exs_pos;
    int curr_i = 0;

    for(std::map<std::string, std::list<std::string> >::iterator it1=genes.begin(); it1!=genes.end(); ++it1) {
        //Forall gene
        for(std::list<std::string>::iterator it2=it1->second.begin(); it2!=it1->second.end(); ++it2) {
            //Forall transcript in gene
            int last_i = -1;
            for(std::list<std::string>::iterator it3=transcripts[*it2].begin(); it3!=transcripts[*it2].end(); ++it3) {
                //Forall exon in transcript
                std::string exon_ID = *it3;
                Feature e = exons[exon_ID];
                std::string pos_ID = getExonID(e.start, e.end);
                try {
                    addedExons.at(pos_ID);
                } catch(const std::out_of_range& oor) {
                    addedExons.insert(std::make_pair(pos_ID, ""));
                    ids_to_index[e.id] = ex_id;
                    std::string curr_ex_string = genomic.substr(e.start-1, e.end-e.start+1);
                    T += curr_ex_string + "|";
                    Exons[ex_id] = curr_ex_string;
                    Exons_Pos.push_back(std::make_pair(e.start, e.end));
                    ++ex_id;
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
    for(const std::pair<int,int>& p1 : Exons_Pos) {
        int j = 1;
        for(const std::pair<int,int>& p2 : Exons_Pos) {
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
    if(edges[x][y] >= 1) {
        return true;
    } else {
        return false;
    }
}

void SplicingGraph::print() const {
    std::cout << T << std::endl;
    for(std::vector<int> v : edges) {
        for(const int& e : v) {
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
    //std::cout << "Saving SG to disk..." << std::endl;
    std::ofstream ofile;
    ofile.open(path + ".sg");
    ofile << reference << " " << ref_length << "\n";
    ofile << T << "\n";
    ofile << exsN << "\n";
    for(const std::vector<int>& v : edges) {
        for(const int& e : v) {
            ofile << e << " ";
        }
        ofile << "; ";
    }
    ofile << "\n";
    for(const std::pair<int,int>& p : Exons_Pos) {
        ofile << p.first << "," << p.second << " ";
    }
    ofile << "\n";
    ofile.close();
}

void SplicingGraph::load(const std::string path) {
    //std::cout << "Loading SG from disk..." << std::endl;
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
