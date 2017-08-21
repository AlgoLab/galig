#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <utility>
#include <list>

#include "FastaReader.hpp"
#include "SplicingGraph.hpp"
#include "bMEM.hpp"
#include "utils.hpp"

#include "MEMsGraph.hpp"

/**
   ### TIME ###
   clock_t t1,t2;
   float t;
   t1=clock();
   ...
   t2=clock();
   t = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
   printf("%.6f,", t);
**/


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

int main(int argc, char* argv[]) {
    std::string genomic;
    std::string annotation;
    std::string rna_seqs;
    int L;
    int eps;
    std::string out1;
    std::string out2;
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
            out1 = std::string(optarg);
            out2 = std::string(optarg) + ".ins";
            break;
        case 'v':
            verbose = true;
            break;
        default:
            printHelp();
            exit(EXIT_FAILURE);
        }
    }
    if(out1.compare("") == 0) {
        out1 = "out.mem";
        out2 = "out.ins";
    }

    // - Building splicing graph
    // --------------------------
    SplicingGraph sg (genomic, annotation);
    if(verbose) sg.print();
    if(index_only) return 0;
    int exsN = sg.getExonsNumber();

    // - Setting up MEMs Index, Reads reader and out file
    // ---------------------------------------------------
    BackwardMEM bm (sg.getText(), genomic);
    FastaReader fastas (rna_seqs);
    std::ofstream outFile1;
    outFile1.open(out1);
    std::ofstream outFile2;
    outFile2.open(out2);

    // - Main loop: one iteration, one read
    // ---------------------------------------
    while(fastas.hasReads()) {
        // Read: (header, read)
        // Path: (weight, [mems])
        std::pair<std::string, std::string> seq = fastas.pop();

        // Original read
        std::string read = seq.second;
        std::list<Mem> mems = bm.getMEMs(read,L);
        std::pair<int, std::list<Mem> > path;
        if(!mems.empty()) {
            MemsGraph mg (read, L, eps, exsN, verbose);
            mg.build(sg, mems);
            path = mg.visit(sg);
        }

        // Reversed-and-complemented read
        std::string readRC = reverseAndComplement(read);
        std::list<Mem> memsRC = bm.getMEMs(readRC,L);
        std::pair<int, std::list<Mem> > pathRC;
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
            outFile1 << strand << " " << seq.first << " " << err << " ";
            for(std::list<Mem>::iterator m=path.second.begin(); m!=path.second.end(); ++m) {
                outFile1 << m->toStr() << " ";
            }
            outFile1 << "\n";
        }
    }
    outFile1.close();
    outFile2.close();
}
