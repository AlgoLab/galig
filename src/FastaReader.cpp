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
                        sequences.push_front(curr_seq);
                        curr_seq = "";
                    }
                    descriptions.push_front(line.substr(1,line.size()-1));
                }
                else {
                    curr_seq += line;
                }
            }
        }
        sequences.push_front(curr_seq);
    } else {
        std::cerr << "Fasta file not found!" << std::endl;
    }
}

int FastaReader::getSize() {
    return descriptions.size();
}

std::pair<std::string, std::string> FastaReader::pop() {
    std::string desc = descriptions.front();
    std::string seq = sequences.front();
    descriptions.pop_front();
    sequences.pop_front();
    return std::make_pair(desc, seq);
}

bool FastaReader::hasReads() {
    return !descriptions.empty();
}
