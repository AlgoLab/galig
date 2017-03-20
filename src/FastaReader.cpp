#include "FastaReader.hpp"

FastaReader::FastaReader(const std::string& fasta) {
    std::ifstream fastaFile (fasta);
    std::string line = "";
    std::string curr_seq = "";
    if(fastaFile.is_open()) {
        while(getline(fastaFile,line)) {
            if(line.compare("") != 0) {
                if(line[0] == '>') {
                    if(curr_seq.compare("") != 0) {
                        sequences.push_back(curr_seq);
                        curr_seq = "";
                    }
                    descriptions.push_back(line.substr(1,line.size()-1));
                }
                else {
                    curr_seq += line;
                }
            }
        }
        sequences.push_back(curr_seq);
    } else {
        std::cerr << "Fasta file not found!" << std::endl;
    }
}

int FastaReader::getSize() {
    return descriptions.size();
}

std::pair<std::string, std::string> FastaReader::getEntry(const int& i) {
    try {
        return std::make_pair(descriptions.at(i), sequences.at(i));
    } catch (const std::out_of_range& e) {
        std::cerr << "Fasta entry exception!" << std::endl;
    }
    return std::make_pair("", "");
}
