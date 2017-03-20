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
    std::cout << "Usage: SGAL [mode] [options]\n" << std::endl;
    std::cout << "Mode:" << std::endl;
    std::cout << "  -1, --index <SG_path>: indexing (it requires -g, -a)" << std::endl;
    std::cout << "  -2, --align <SG_path>: aligning (it requires -r, -l, -k)" << std::endl << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -g, --genomic <path>:" << std::endl;
    std::cout << "  -a, --annotation <path>" << std::endl;
    std::cout << "  -r, --reads <path>" << std::endl;
    std::cout << "  -l, --L <int>: MEMs length" << std::endl;
    std::cout << "  -e, --eps <int>: " << std::endl;
    std::cout << "  -v, --verbose:" << std::endl;
}

int main(int argc, char* argv[]) {
    bool mode;
    std::string sg_index;
    std::string genomic;
    std::string annotation;
    std::string rna_seqs;
    int L;
    int eps;

    int c;
    while (1) {
        static struct option long_options[] =
            {
                {"index", no_argument, 0, '1'},
                {"align", no_argument, 0, '2'},
                {"verbose", no_argument, 0, 'v'},
                {"genomic", required_argument, 0, 'g'},
                {"annotation", required_argument, 0, 'a'},
                {"reads",  required_argument, 0, 'r'},
                {"L",  required_argument, 0, 'l'},
                {"eps",    required_argument, 0, 'e'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "1:2:g:a:r:l:e:v", long_options, &option_index);

        if (c == -1) {
            break;
        }

        switch(c) {
        case '1':
            mode = true;
            sg_index = optarg;
            break;
        case '2':
            mode = false;
            sg_index = optarg;
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
        case 'v':
            std::cout << "verbose" << std::endl;
            break;
        default:
            printHelp();
            exit(EXIT_FAILURE);
        }
    }
    //if(mode) {
    SplicingGraph sg (genomic, annotation, sg_index);
    sg.print();
    //exit(0);
    //}
    //else {
    //SplicingGraph sg (sg_index);
    FastaReader fastas (rna_seqs);
    BackwardMEM bm (sg.getText(), genomic);
    int i = 0;

    std::ofstream outFile;
    outFile.open("../OUT.mem");
    while(i<fastas.getSize()) {
        std::pair<std::string, std::string> seq = fastas.getEntry(i);
        std::string read = seq.second;
        std::list<Mem> mems = bm.getMEMs(read,L);
        // std::cout << seq.first << " ";
        // for(Mem m : mems) {
        //     std::cout << m.toStr() << " ";
        // }
        // std::cout << std::endl;
        MemsGraph mg (sg, read, mems, L, eps);
        mg.build(sg, mems);
        std::pair<int, std::list<std::list<Mem> > > paths = mg.visit();

        for(std::list<std::list<Mem> >::iterator p=paths.second.begin(); p!=paths.second.end(); ++p) {
          outFile << "+ " << seq.first << " " << paths.first << " ";
          for(std::list<Mem>::iterator m=p->begin(); m!=p->end(); ++m) {
            outFile << m->toStr() << " ";
          }
          outFile << "\n";
        }

        /**
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
        std::cout << flag << std::endl;
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
        **/
        ++i;
        if(i%100 == 0) {
            std::cout << "Processed " << i << " reads." << std::endl;
        }
    }
    outFile.close();
}
