#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>

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
    std::cout << "  -k, --K <int>: " << std::endl;
    std::cout << "  -v, --verbose:" << std::endl;
}

int main(int argc, char* argv[]) {
    bool mode;
    std::string sg_index;
    std::string genomic;
    std::string annotation;
    std::string rna_seqs;
    int L;
    int K;

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
                {"K",    required_argument, 0, 'k'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "1:2:g:a:r:l:k:v", long_options, &option_index);

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
        case 'k':
            K = std::stoi(optarg);
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
    //sg.print();
    //exit(0);
    //}
    //else {
    //SplicingGraph sg (sg_index);
    //sg.print();
    FastaReader fastas (rna_seqs);
    BackwardMEM bm (sg.getText(), genomic);
    int i = 0;
    //std::ofstream outFile;
    //outFile.open(genomic + "_res");
    while(i<fastas.getSize()) {
        std::pair<std::string, std::string> seq = fastas.getEntry(i);
        std::string read = seq.second;
        std::list<Mem> mems = bm.getMEMs(read,L);
        MemsGraph mg (sg, mems, L);
        ++i;
    }
}
        /**
           MemsList ml (seq.second.size());
           for(const Mem& mem : mems) {
           ml.addMem(mem);
           }
           MemsGraph mg (sg, ml, K, 80);
           //mg.saveImage(seq.first);
           //std::cout << "." << std::endl;
           //mg.print();
           mg.visit();
           std::list<std::string> out (mg.getOutput());
           ++i;
           //for(std::string s : out) {
           //    outFile << seq.first << " " << s << "\n";
           //}
           if(i%50 == 0) {
           std::cout << "Processed " << i << " reads." << std::endl;
           }
        **/
    //outFile.close();
