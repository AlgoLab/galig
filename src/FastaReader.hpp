#ifndef _FASTA_READER_HPP_
#define _FASTA_READER_HPP_

#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <utility>
#include <iterator>

class FastaReader {
private:
    std::list<std::string> descriptions;
    std::list<std::string> sequences;
public:
    FastaReader(const std::string&);
    int getSize();
    std::pair<std::string, std::string> pop();
};
#endif
