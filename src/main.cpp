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

void printHelp() {
    std::cout << "Usage: SGAL [options] (required: -g -a -r -l -e)\n" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -g, --genomic <path>:" << std::endl;
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
    std::string out;
    bool verbose = false;

    int c;
    while (1) {
        static struct option long_options[] =
            {
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
        c = getopt_long(argc, argv, "g:a:r:l:e:o:v", long_options, &option_index);

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
            out = optarg;
            break;
        case 'v':
            verbose = true;
            break;
        default:
            printHelp();
            exit(EXIT_FAILURE);
        }
    }
    if(verbose) {
        std::cout << "Starting..." << std::endl;
    }
    SplicingGraph sg (genomic, annotation);
    if(verbose) {
        sg.print();
    }
    BackwardMEM bm (sg.getText(), genomic);
    FastaReader fastas (rna_seqs);
    int i = 0;

    if(out.compare("") == 0) {
        out = "OUT.mem";
    }
    std::ofstream outFile;
    outFile.open(out);
    while(i<fastas.getSize()) {
        std::pair<std::string, std::string> seq = fastas.getEntry(i);

        std::string read = seq.second;
        std::pair<int, std::list<std::list<Mem> > > paths;
        std::list<Mem> mems = bm.getMEMs(read,L);
        bool flag = false;
        if(mems.size() != 0) {
            MemsGraph mg (sg, read, L, eps, verbose);
            mg.build(sg, mems);
            paths = mg.visit();
            flag = true;
        }

        std::string read_RC = reverse_and_complement(read);
        std::pair<int, std::list<std::list<Mem> > > paths_RC;
        std::list<Mem> mems_RC = bm.getMEMs(read_RC,L);
        bool flag_RC = false;
        if(mems_RC.size() != 0) {
            MemsGraph mg_RC (sg, read_RC, L, eps, verbose);
            mg_RC.build(sg, mems_RC);
            paths_RC = mg_RC.visit();
            flag_RC = true;
        }
        int best = 0;
        if(flag && !flag_RC) {
            best = 1;
        } else if(!flag && flag_RC) {
            best = 2;
        } else {
            if(paths.first <= paths_RC.first) {
                best = 1;
            } else {
                best = 2;
            }
        }
        switch(best) {
        case 0:
            break;
        case 1:
            for(std::list<std::list<Mem> >::iterator p=paths.second.begin(); p!=paths.second.end(); ++p) {
                outFile << "+ " << seq.first << " " << paths.first << " ";
                for(std::list<Mem>::iterator m=p->begin(); m!=p->end(); ++m) {
                    outFile << m->toStr() << " ";
                }
                outFile << "\n";
            }
            break;
        case 2:
            for(std::list<std::list<Mem> >::iterator p=paths_RC.second.begin(); p!=paths_RC.second.end(); ++p) {
                outFile << "- " << seq.first << " " << paths_RC.first << " ";
                for(std::list<Mem>::iterator m=p->begin(); m!=p->end(); ++m) {
                    outFile << m->toStr() << " ";
                }
                outFile << "\n";
            }
            break;
        }
        ++i;
        if(i%1000 == 0) {
            std::cout << "Processed " << i << " reads." << std::endl;
        }
    }
    outFile.close();
    if(verbose) {
        std::cout << "...Ending" << std::endl;
    }
}
