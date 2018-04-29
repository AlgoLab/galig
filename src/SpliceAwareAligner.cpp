#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <utility>
#include <list>

#include <zlib.h>
#include <stdio.h>

#include "kseq.h"

#include "SplicingGraph.hpp"
#include "bMEM.hpp"
#include "utils.hpp"
#include "MEMsGraph.hpp"

KSEQ_INIT(gzFile, gzread)

void printHelp() {
    std::cout << "Usage: SGAL [options] (required: -g -a -s -o)\n" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -g, --genome <path>" << std::endl;
    std::cout << "  -a, --annotation <path>" << std::endl;
    std::cout << "  -s, --sample <path>" << std::endl;
    std::cout << "  -o, --output <path>: output file" << std::endl;
    std::cout << "  -l, --L <int>: minimum lenght of MEMs used to build the alignments (default: 15)" << std::endl;
    std::cout << "  -e, --eps <int>: error rate, a value from 0 to 100 (default: 3)" << std::endl;
    std::cout << "  -h, --help: show this help message and exit" << std::endl;
    //std::cout << "  -v, --verbose: explain what is being done and save .dot" << std::endl;
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
    if(verbose) {
        for(const Mem& m : mems) {
            std::cout << m.toStr() << " ";
        }
        std::cout << std::endl;
    }
    std::list<std::pair<int, std::list<Mem> > > paths; // Path: [(weight, [mems])]
    if(!mems.empty()) {
        MemsGraph mg (read, L, eps, exsN, verbose);
        mg.build(sg, mems);
        paths = mg.visit(sg);
    }
    // Reversed-and-complemented read
    std::string readRC = reverseAndComplement(read);
    std::list<Mem> memsRC = bm.getMEMs(readRC,L);
    if(verbose) {
        for(const Mem& m : memsRC) {
            std::cout << m.toStr() << " ";
        }
        std::cout << std::endl;
    }
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
    int L = 0;
    int eps = -1;
    std::string out;
    bool verbose = false;

    // - Collecting command line parameters
    // -------------------------------------
    int c;
    while (1) {
        static struct option long_options[] =
            {
                {"genomic", required_argument, 0, 'g'},
                {"annotation", required_argument, 0, 'a'},
                {"sample",  required_argument, 0, 's'},
                {"L",  required_argument, 0, 'l'},
                {"erate",    required_argument, 0, 'e'},
                {"output", required_argument, 0, 'o'},
                {"help", no_argument, 0, 'h'},
                //{"verbose", no_argument, 0, 'v'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "g:a:s:l:e:o:h", long_options, &option_index);

        if (c == -1) {
            break;
        }

        switch(c) {
        case 'g':
            genomic = optarg;
            break;
        case 'a':
            annotation = optarg;
            break;
        case 's':
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
        // case 'v':
        //     verbose = true;
        //     break;
        case 'h':
            printHelp();
            exit(EXIT_SUCCESS);
        default:
            printHelp();
            exit(EXIT_FAILURE);
        }
    }
    if(L == 0) {
        L = 15;
    }
    if(eps < 0 || eps > 100) {
        eps = 3;
    }

    // - Building splicing graph
    // --------------------------
    gzFile fastain = gzopen(genomic.c_str(), "r");
    kseq_t *reference = kseq_init(fastain);
    kseq_read(reference);

    SplicingGraph sg (reference->seq.s, annotation);

    // sg.print();
    kseq_destroy(reference);
    gzclose(fastain);

    int exsN = sg.getExonsNumber();

    // - Setting up MEMs Index and out file
    // ---------------------------------------------------
    BackwardMEM bm (sg.getText(), genomic);

    std::ofstream outFile;
    outFile.open(out);

    fastain = gzopen(rna_seqs.c_str(), "r");
    std::pair<char, std::list<std::pair<int, std::list<Mem> > > > paths;

    kseq_t *seqs = kseq_init(fastain);
    int l;
    std::string head;
    std::string read;
    int i = 1;

    // - Main loop: one iteration, one read
    // ---------------------------------------
    while ((l = kseq_read(seqs)) >= 0) {
        head = seqs->name.s;
        read = seqs->seq.s;

        paths = analyzeRead(bm, sg, read, L, eps, exsN, verbose);
        if(paths.first != '/') {
            for(std::pair<int, std::list<Mem> > path : paths.second) {
                if(!path.second.empty()) {
                    int err = path.first;
                    outFile << paths.first << " " << head << " " << err << " ";
                    for(std::list<Mem>::iterator m=path.second.begin(); m!=path.second.end(); ++m) {
                        outFile << m->toStr() << " ";
                    }
                    if(paths.first == '+')
                        outFile << read;
                    else
                        outFile << reverseAndComplement(read);
                    outFile << "\n";
                }
            }
        }
        if(i%10000 == 0)
            std::cout << "Processed " << i << " genes." << std::endl;
        ++i;
    }
    kseq_destroy(seqs);
    gzclose(fastain);
    return 0;
}
