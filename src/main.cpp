#include <iostream>
#include <string>
#include <fstream>

#include "FastaReader.hpp"
#include "SplicingGraph.hpp"
#include "bMEM.hpp"
#include "utils.hpp"

#include "MEMsList.hpp"
#include "MEMsGraph.hpp"

int main(int argc, char* argv[]) {
    std::string genomic = argv[1];
    std::string annotation = argv[2];
    std::string rna_seqs = argv[3];
    int L = std::stoi(argv[4]);
    int K = std::stoi(argv[5]);

    FastaReader fastas (rna_seqs);
    SplicingGraph sg (genomic, annotation);
    //sg.print();
    BackwardMEM bm (sg.getText(), genomic);
    int i = 0;
    std::ofstream outFile;
    outFile.open(genomic + "_res");
    while(i<fastas.getSize()) {
        std::pair<std::string, std::string> seq = fastas.getEntry(i);
        //std::string x = reverse_and_complement(seq.second);
        std::list<Mem> mems = bm.getMEMs(seq.second,L);
        //bm.printMEMs();
        MemsList ml (seq.second.size());
        for(const Mem& mem : mems) {
            ml.addMem(mem);
        }
        MemsGraph mg (sg, ml, K, 80);
        //mg.saveImage(seq.first);
        mg.visit();
        std::list<std::string> out (mg.getOutput());
        ++i;

        for(std::string s : out) {
            outFile << seq.first << " " << s << "\n";
        }
    }
    outFile.close();
}
