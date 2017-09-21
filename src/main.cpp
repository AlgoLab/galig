#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <utility>
#include <list>

#include "SplicingGraph.hpp"
#include "bMEM.hpp"
#include "utils.hpp"

#include "MEMsGraph.hpp"

void printHelp() {
    std::cout << "Usage: SGAL [options] (required: -g -a -r -l -e)\n" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -I, --index: index only" << std::endl;
    std::cout << "  -g, --genomic <path>" << std::endl;
    std::cout << "  -a, --annotation <path>" << std::endl;
    std::cout << "  -r, --reads <path>" << std::endl;
    std::cout << "  -l, --L <int>: MEMs length" << std::endl;
    std::cout << "  -e, --eps <int>: " << std::endl;
    std::cout << "  -o, --output <path>: output path" << std::endl;
    std::cout << "  -v, --verbose: explain what is being done and save .dot" << std::endl;
}

std::pair<char, std::list<std::pair<int, std::list<Mem> > > > analyzeRead(BackwardMEM& bm,
                                                                                const SplicingGraph& sg,
                                                                                const std::string& read,
                                                                                const int& L,
                                                                                const int& eps,
                                                                                const int& exsN,
                                                                                const bool& verbose) {
    // Original read
    std::list<Mem> mems = bm.getMEMs(read,L);
    std::list<std::pair<int, std::list<Mem> > > paths; // Path: [(weight, [mems])]
    if(!mems.empty()) {
        MemsGraph mg (read, L, eps, exsN, verbose);
        mg.build(sg, mems);
        paths = mg.visit(sg);
    }
    // Reversed-and-complemented read
    std::string readRC = reverseAndComplement(read);
    std::list<Mem> memsRC = bm.getMEMs(readRC,L);
    std::list<std::pair<int, std::list<Mem> > > pathsRC; // Path: [(weight, [mems])]
    if(!memsRC.empty()) {
        MemsGraph mgRC (readRC, L, eps, exsN, verbose);
        mgRC.build(sg, memsRC);
        pathsRC = mgRC.visit(sg);
    }
    bool empty = paths.empty();
    bool emptyRC = pathsRC.empty();
    char strand = '/';
    if(!empty || !emptyRC) {
        if(!empty && emptyRC) {
            strand = '+';
        } else if(empty && !emptyRC) {
            paths = pathsRC;
            strand = '-';
        } else {
            if(paths.front().first <= pathsRC.front().first) {
                strand = '+';
            } else {
                paths = pathsRC;
                strand = '-';
            }
        }
    }
    return std::make_pair(strand, paths);
}

int main(int argc, char* argv[]) {
    std::string genomic;
    std::string annotation;
    std::string rna_seqs;
    int L;
    int eps;
    std::string out;
    bool index_only = false;
    bool verbose = false;

    // - Collecting command line parameters
    // -------------------------------------
    int c;
    while (1) {
        static struct option long_options[] =
            {
                {"index", no_argument, 0, 'I'},
                {"genomic", required_argument, 0, 'g'},
                {"annotation", required_argument, 0, 'a'},
                {"reads",  required_argument, 0, 'r'},
                {"L",  required_argument, 0, 'l'},
                {"eps",    required_argument, 0, 'e'},
                {"output", required_argument, 0, 'o'},
                {"verbose", no_argument, 0, 'v'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "Ig:a:r:l:e:o:v", long_options, &option_index);

        if (c == -1) {
            break;
        }

        switch(c) {
        case 'I':
            index_only = true;
            break;
        case 'g':
            genomic = optarg;
            break;
        case 'a':
            annotation = optarg;
            break;
        case 'r':
            rna_seqs = optarg;
            break;
        case 'l':
            L = std::stoi(optarg);
            break;
        case 'e':
            eps = std::stoi(optarg);
            break;
        case 'o':
            out = std::string(optarg);
            break;
        case 'v':
            verbose = true;
            break;
        default:
            printHelp();
            exit(EXIT_FAILURE);
        }
    }
    if(out.compare("") == 0) {
        out = "out.mem";
    }

    // - Building splicing graph
    // --------------------------
    SplicingGraph sg (genomic, annotation);
    if(verbose) sg.print();
    if(index_only) return 0;
    int exsN = sg.getExonsNumber();

    // - Setting up MEMs Index and out file
    // ---------------------------------------------------
    BackwardMEM bm (sg.getText(), genomic);
    std::ofstream outFile;
    outFile.open(out);

    std::ifstream fastaFile (rna_seqs);
    std::pair<char, std::list<std::pair<int, std::list<Mem> > > > paths;
    // - Main loop: one iteration, one read
    // ---------------------------------------
    if(fastaFile.is_open()) {
        std::string line;
        std::string head;
        std::string read;
        while(getline(fastaFile,line)) {
            if(line.compare("") != 0) {
                if(line[0] == '>') {
                    if(read.compare("") != 0) {
                        paths = analyzeRead(bm, sg, read, L, eps, exsN, verbose);
                        if(paths.first != '/') {
                            for(std::pair<int, std::list<Mem> > path : paths.second) {
                                if(!path.second.empty()) {
                                    int err = path.first;
                                    outFile << paths.first << " " << head << " " << err << " ";
                                    for(std::list<Mem>::iterator m=path.second.begin(); m!=path.second.end(); ++m) {
                                        outFile << m->toStr() << " ";
                                    }
                                    outFile << "\n";
                                }
                            }
                        }
                        read = "";
                    }
                    head = line.substr(1,line.size()-1);
                } else {
                    read += line;
                }
            }
        }
        paths = analyzeRead(bm, sg, read, L, eps, exsN, verbose);
        if(paths.first != '/') {
            for(std::pair<int, std::list<Mem> > path : paths.second) {
                if(!path.second.empty()) {
                    int err = path.first;
                    outFile << paths.first << " " << head << " " << err << " ";
                    for(std::list<Mem>::iterator m=path.second.begin(); m!=path.second.end(); ++m) {
                        outFile << m->toStr() << " ";
                    }
                    outFile << "\n";
                }
            }
        }
    } else {
        std::cerr << "Reads file not found!" << std::endl;
    }
    outFile.close();
}
