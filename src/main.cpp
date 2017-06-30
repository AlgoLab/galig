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
    std::cout << "  -g, --genomic <path>:" << std::endl;
    std::cout << "  -a, --annotation <path>" << std::endl;
    std::cout << "  -r, --reads <path>" << std::endl;
    std::cout << "  -l, --L <int>: MEMs length" << std::endl;
    std::cout << "  -e, --eps <int>: " << std::endl;
    std::cout << "  -o, --output <path>: output path" << std::endl;
    std::cout << "  -G, --greedy: run greedy search instead of exhaustive one" << std::endl;
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
    bool greedy = true;

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
                {"exhaustive", no_argument, 0, 'X'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "g:a:r:l:e:o:vX", long_options, &option_index);

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
        case 'X':
            greedy = false;
            break;
        default:
            printHelp();
            exit(EXIT_FAILURE);
        }
    }

    SplicingGraph sg (genomic, annotation);
    int exsN = sg.getExonsNumber();
    if(verbose) {
        sg.print();
    }
    BackwardMEM bm (sg.getText(), genomic);
    FastaReader fastas (rna_seqs);

    int i = 0;
    if(out.compare("") == 0) {
        out = "out.mem";
    }
    std::ofstream outFile;
    outFile.open(out);
    while(i<fastas.getSize()) {
        std::pair<std::string, std::string> seq = fastas.getEntry(i);
        std::string read = seq.second;
        std::string read_RC = reverse_and_complement(read);
        std::list<Mem> mems = bm.getMEMs(read,L);
        std::list<Mem> mems_RC = bm.getMEMs(read_RC,L);

        if(greedy) {
            bool flag = false;
            bool flag_RC = false;
            std::pair<int, std::list<Mem> > res;
            std::pair<int, std::list<Mem> > res_RC;

            if(!mems.empty()) {
                if(verbose) {
                    for(const Mem& m : mems) {
                        std::cout << m.toStr() << " " << std::endl;
                    }
                    std::cout << std::endl;
                }
                MemsGraph mg (read, L, eps, exsN, verbose);
                res = mg.build_greedy(sg, mems);
                flag = true;
            }
            if(!mems_RC.empty()) {
                if(verbose) {
                    for(const Mem& m : mems) {
                        std::cout << m.toStr() << " " << std::endl;
                    }
                    std::cout << std::endl;
                }
                MemsGraph mg_RC (read_RC, L, eps, exsN, verbose);
                res_RC = mg_RC.build_greedy(sg, mems_RC);
                flag_RC = true;
            }
            int best = 0;
            if(flag && !flag_RC) {
                if(!res.second.empty()) {
                    best = 1;
                }
            } else if(!flag && flag_RC) {
                if(!res_RC.second.empty()) {
                    best = 2;
                }
            } else {
                if(res.first < res_RC.first) {
                    if(!res.second.empty()) {
                        best = 1;
                    }
                } else if(res.first > res_RC.first) {
                    if(!res_RC.second.empty()) {
                        best = 2;
                    }
                } else {
                    if(!res.second.empty() && res_RC.second.empty()) {
                        best = 1;
                    } else if(res.second.empty() && !res_RC.second.empty()) {
                        best = 2;
                    } else {
                        
                    }
                }
            }
            switch(best) {
            case 0:
                break;
            case 1:
                outFile << "+ " << seq.first << " " << res.first << " ";
                for(std::list<Mem>::iterator m=res.second.begin(); m!=res.second.end(); ++m) {
                    outFile << m->toStr() << " ";
                }
                outFile << "\n";
                break;
            case 2:
                outFile << "- " << seq.first << " " << res_RC.first << " ";
                for(std::list<Mem>::iterator m=res_RC.second.begin(); m!=res_RC.second.end(); ++m) {
                    outFile << m->toStr() << " ";
                }
                outFile << "\n";
                break;
            }
            
        } else {
            bool flag = false;
            bool flag_RC = false;
            std::pair<bool, std::pair<int, std::list<std::pair<bool, std::list<Mem> > > > > paths;
            std::pair<bool, std::pair<int, std::list<std::pair<bool, std::list<Mem> > > > > paths_RC;

            if(!mems.empty()) {
                MemsGraph mg (read, L, eps, exsN, verbose);
                mg.build(sg, mems);
                paths = mg.visit(sg);
                flag = true;
            }
            if(!mems_RC.empty()) {
                MemsGraph mg_RC (read_RC, L, eps, exsN, verbose);
                mg_RC.build(sg, mems_RC);
                paths_RC = mg_RC.visit(sg);
                flag_RC = true;
            }

            int best = 0;
            if(flag && !flag_RC) {
                best = 1;
            } else if(!flag && flag_RC) {
                best = 2;
            } else {
                if(paths.second.first <= paths_RC.second.first) {
                    best = 1;
                } else {
                    best = 2;
                }
            }

            bool atLeastOneAnnotaded;
            int err;
            switch(best) {
            case 0:
                break;
            case 1:
                atLeastOneAnnotaded = paths.first;
                err = paths.second.first;
                for(std::list<std::pair<bool, std::list<Mem> > >::iterator p=paths.second.second.begin(); p!=paths.second.second.end(); ++p) {
                    bool annotated = p->first;
                    bool f = true;
                    if(atLeastOneAnnotaded and !annotated) {
                        f = false;
                    }
                    if(f) {
                        outFile << "+ " << seq.first << " " << err << " ";
                        for(std::list<Mem>::iterator m=p->second.begin(); m!=p->second.end(); ++m) {
                            outFile << m->toStr() << " ";
                        }
                        outFile << "\n";
                    }
                }
                break;
            case 2:
                atLeastOneAnnotaded = paths_RC.first;
                err = paths_RC.second.first;
                for(std::list<std::pair<bool, std::list<Mem> > >::iterator p=paths_RC.second.second.begin(); p!=paths_RC.second.second.end(); ++p) {
                    bool annotated = p->first;
                    bool f = true;
                    if(atLeastOneAnnotaded and !annotated) {
                        f = false;
                    }
                    if(f && p->second.size()>0) {
                        outFile << "- " << seq.first << " " << err << " ";
                        for(std::list<Mem>::iterator m=p->second.begin(); m!=p->second.end(); ++m) {
                            outFile << m->toStr() << " ";
                        }
                        outFile << "\n";
                    }
                }
                break;
            }
        }
        ++i;
        // if(i%1000 == 0) {
        //     std::cout << "Processed " << i << " reads." << std::endl;
        // }
    }
    outFile.close();
    return 0;
}
