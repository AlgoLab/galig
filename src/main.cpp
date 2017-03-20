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
}

int main(int argc, char* argv[]) {
    std::string sg_index;
    std::string genomic;
    std::string annotation;
    std::string rna_seqs;
    int L;
    int eps;
    std::string out = "";

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
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "g:a:r:l:e:o", long_options, &option_index);

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
            out = std::stoi(optarg);
            break;
        default:
            printHelp();
            exit(EXIT_FAILURE);
        }
    }

    SplicingGraph sg (genomic, annotation);
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
        
        std::list<Mem> mems = bm.getMEMs(read,L);

        MemsGraph mg (sg, read, mems, L, eps);
        mg.build(sg, mems);
        std::pair<int, std::list<std::list<Mem> > > paths = mg.visit();

        std::string read1 = reverse_and_complement(read);
        std::list<Mem> mems1 = bm.getMEMs(read1,L);
        MemsGraph mg1 (sg, read1, mems, L, eps);
        mg1.build(sg, mems1);
        std::pair<int, std::list<std::list<Mem> > > paths1 = mg1.visit();

        int flag = 0;
        if(paths.first <= paths1.first) {
            if(paths.second.size() != 0) {
                flag = 1;
            } else if(paths1.second.size() != 0) {
                flag = 2;
            }
        } else {
            if(paths1.second.size() != 0) {
                flag = 2;
            } else if(paths1.second.size() != 0) {
                flag = 1;
            }
        }

        switch(flag) {
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
            for(std::list<std::list<Mem> >::iterator p=paths1.second.begin(); p!=paths1.second.end(); ++p) {
                outFile << "- " << seq.first << " " << paths.first << " ";
                for(std::list<Mem>::iterator m=p->begin(); m!=p->end(); ++m) {
                    outFile << m->toStr() << " ";
                }
                outFile << "\n";
            }
            break;
        }
        ++i;
        if(i%100 == 0) {
            std::cout << "Processed " << i << " reads." << std::endl;
        }
    }
    outFile.close();
}
