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
    bool twoPass = false;
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
                {"2pass", required_argument, 0, 't'},
                {"verbose", no_argument, 0, 'v'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "Ig:a:r:l:e:o:tv", long_options, &option_index);

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
            out1 = std::string(optarg) + ".mem";
            out2 = std::string(optarg) + ".ins";
            break;
        case 't':
            twoPass = true;
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
    int i = 0;
    int numReads = fastas.getSize();
    while(i<numReads) {
        // Read: (header, read)
        // Paths: <<atLeastOneAnnotated, weight>, [<annotated, [mems]>]>
        std::pair<std::string, std::string> seq = fastas.pop();

        // Original read
        std::string read = seq.second;
        std::list<Mem> mems = bm.getMEMs(read,L);
        std::pair<std::pair<bool, int>, std::list<std::pair<bool, std::list<Mem> > > > paths;
        if(!mems.empty()) {
            MemsGraph mg (read, L, eps, exsN, verbose);
            mg.build(sg, mems, false);
            paths = mg.visit(sg);
        }

        // Reversed-and-complemented read
        std::string read_RC = reverse_and_complement(read);
        std::list<Mem> mems_RC = bm.getMEMs(read_RC,L);
        std::pair<std::pair<bool, int>, std::list<std::pair<bool, std::list<Mem> > > > paths_RC;
        if(!mems_RC.empty()) {
            MemsGraph mg_RC (read_RC, L, eps, exsN, verbose);
            mg_RC.build(sg, mems_RC, false);
            paths_RC = mg_RC.visit(sg);
        }

        bool empty = paths.second.empty();
        bool empty_RC = paths_RC.second.empty();
        char strand;
        if(empty && empty_RC) {
            if(twoPass) {
                MemsGraph mg (read, L, eps, exsN, verbose);
                mg.build(sg, mems, true);
                paths = mg.visit(sg);
                MemsGraph mg_RC (read_RC, L, eps, exsN, verbose);
                mg_RC.build(sg, mems_RC, true);
                paths_RC = mg_RC.visit(sg);
                for(std::list<std::pair<bool, std::list<Mem> > >::iterator p=paths.second.begin(); p!=paths.second.end(); ++p) {
                    outFile2 << "+ " << seq.first << " 0 ";
                    for(std::list<Mem>::iterator m=p->second.begin(); m!=p->second.end(); ++m) {
                        outFile2 << m->toStr() << " ";
                    }
                    outFile2 << "\n";
                }
                for(std::list<std::pair<bool, std::list<Mem> > >::iterator p=paths_RC.second.begin(); p!=paths_RC.second.end(); ++p) {
                    outFile2 << "- " << seq.first << " 0 ";
                    for(std::list<Mem>::iterator m=p->second.begin(); m!=p->second.end(); ++m) {
                        outFile2 << m->toStr() << " ";
                    }
                    outFile2 << "\n";
                }
            }
        } else {
            if(!empty && empty_RC) {
                strand = '+';
            } else if(empty && !empty_RC) {
                paths = paths_RC;
                strand = '-';
            } else {
                if(paths.first.second <= paths_RC.first.second) {
                    strand = '+';
                } else {
                    paths = paths_RC;
                    strand = '-';
                }
            }

            bool atLeastOneAnnotaded = paths.first.first;
            int err = paths.first.second;
            for(std::list<std::pair<bool, std::list<Mem> > >::iterator p=paths.second.begin(); p!=paths.second.end(); ++p) {
                bool annotated = p->first;
                bool f = true;
                if(atLeastOneAnnotaded and !annotated) {
                    f = false;
                }
                if(f) {
                    outFile1 << strand << " " << seq.first << " " << err << " ";
                    for(std::list<Mem>::iterator m=p->second.begin(); m!=p->second.end(); ++m) {
                        outFile1 << m->toStr() << " ";
                    }
                    outFile1 << "\n";
                    if(annotated)
                        break;
                }
            }
        }
        ++i;
        // if(i%1000 == 0) {
        //     std::cout << "Processed " << i << " reads." << std::endl;
        // }
    }
    outFile1.close();
    outFile2.close();
}
