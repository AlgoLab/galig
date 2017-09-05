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

void analyzeRead(BackwardMEM& bm,
                 const SplicingGraph& sg,
                 const std::string& head,
                 const std::string& read,
                 const int& L,
                 const int& eps,
                 const int& exsN,
                 const bool& verbose,
                 std::ofstream& outFile) {
    // Original read
    std::list<Mem> mems = bm.getMEMs(read,L);
    std::pair<int, std::list<Mem> > path; // Path: (weight, [mems])
    if(!mems.empty()) {
        MemsGraph mg (read, L, eps, exsN, verbose);
        mg.build(sg, mems);
        path = mg.visit(sg);
    }
    
    // Reversed-and-complemented read
    std::string readRC = reverseAndComplement(read);
    std::list<Mem> memsRC = bm.getMEMs(readRC,L);
    std::pair<int, std::list<Mem> > pathRC; // Path: (weight, [mems])
    if(!memsRC.empty()) {
        MemsGraph mgRC (readRC, L, eps, exsN, verbose);
        mgRC.build(sg, memsRC);
        pathRC = mgRC.visit(sg);
    }
    int err = path.first;
    int errRC = pathRC.first;
    bool empty = path.second.empty();
    bool emptyRC = pathRC.second.empty();
    char strand;
    if(!empty || !emptyRC) {
        if(!empty && emptyRC) {
            strand = '+';
        } else if(empty && !emptyRC) {
            path = pathRC;
            err = errRC;
            strand = '-';
        } else {
            if(path.first <= pathRC.first) {
                strand = '+';
            } else {
                path = pathRC;
                err = errRC;
                strand = '-';
            }
        }
        outFile << strand << " " << head << " " << err << " ";
        for(std::list<Mem>::iterator m=path.second.begin(); m!=path.second.end(); ++m) {
            outFile << m->toStr() << " ";
        }
        outFile << "\n";
    }
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

    // - Main loop: one iteration, one read
    // ---------------------------------------
    std::ifstream fastaFile (rna_seqs);
    if(fastaFile.is_open()) {
        std::string line;
        std::string head;
        std::string read;
        while(getline(fastaFile,line)) {
            if(line.compare("") != 0) {
                if(line[0] == '>') {
                    if(read.compare("") != 0) {
                        analyzeRead(bm, sg, head, read, L, eps, exsN, verbose, outFile);
                        read = "";
                    }
                    head = line.substr(1,line.size()-1);
                } else {
                    read += line;
                }
            }
        }
        analyzeRead(bm, sg, head, read, L, eps, exsN, verbose, outFile);
    } else {
        std::cerr << "Reads file not found!" << std::endl;
    }
    outFile.close();
}
