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
    std::cout << "Usage: SGAL [options] (required: -g -a -1 -2 -o)\n" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -g, --genome <path>" << std::endl;
    std::cout << "  -a, --annotation <path>" << std::endl;
    std::cout << "  -1, --sample 1 <path>" << std::endl;
    std::cout << "  -2, --sample 2 <path>" << std::endl;
    std::cout << "  -o, --output <path>: output file" << std::endl;
    std::cout << "  -f, --flt : fragment library type (default: ISF)" << std::endl;
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

// Only aligns the provided read, NOT the reverse-complemented version
std::pair<char, std::list<std::pair<int, std::list<Mem> > > > analyzeReadSimple(BackwardMEM& bm,
                                                                                const SplicingGraph& sg,
                                                                                const std::string& read,
                                                                                const int& L,
                                                                                const int& eps,
                                                                                const int& exsN,
                                                                                const bool& verbose,
                                                                                const char& strandIn) {
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

    bool empty = paths.empty();
    char strand = '/';
    if(!empty)
        strand = strandIn;
    return std::make_pair(strand, paths);
}

int main(int argc, char* argv[]) {

    std::string genomic;
    std::string annotation;
    std::string rna_seq_1, rna_seq_2;
    std::string flt="ISF";      //fragment library type
    int L = 0;
    int eps = -1;
    std::string out, out_1, out_2;
    bool verbose = false;

    // - Collecting command line parameters
    // -------------------------------------
    int c;
    while (1) {
        static struct option long_options[] =
            {
                {"genomic", required_argument, 0, 'g'},
                {"annotation", required_argument, 0, 'a'},
                {"sample1",  required_argument, 0, '1'},
                {"sample2",  required_argument, 0, '2'},
                {"L",  required_argument, 0, 'l'},
                {"erate",    required_argument, 0, 'e'},
                {"output", required_argument, 0, 'o'},
                {"flt", required_argument, 0, 'f'},
                {"help", no_argument, 0, 'h'},
                //{"verbose", no_argument, 0, 'v'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "g:a:1:2:l:e:o:f:h", long_options, &option_index);

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
        case '1':
            rna_seq_1 = optarg;
            break;
        case '2':
            rna_seq_2 = optarg;
            break;
        case 'l':
            L = std::stoi(optarg);
            break;
        case 'e':
            eps = std::stoi(optarg);
            break;
        case 'f':
            flt = optarg;
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

    // - Extract fragment library type from options
    bool inward=false, outward=false, matching=false;
    bool stranded=false, unstranded=false;
    bool forward=false, reverse=false;
    for(char &c : flt) {
        switch(c) {
            case 'I':
                inward = true;
                break;
            case 'O':
                outward = true;
                break;
            case 'M':
                matching = true;
                break;
            case 'S':
                stranded = true;
                break;
            case 'U':
                unstranded = true;
                break;
            case 'F':
                forward = true;
                break;
            case 'R':
                reverse = true;
                break;
        }
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

    std::ofstream outFile_1, outFile_2;

    // Suppose out has no extension
    out_1 = out + "-1" + ".mem";
    out_2 = out + "-2" + ".mem";

    outFile_1.open(out_1);
    outFile_2.open(out_2);

    //Declare variables
    std::pair<char, std::list<std::pair<int, std::list<Mem> > > > paths;
    kseq_t *seqs_1, *seqs_2;
    int l1, l2;
    std::string head_1, head_2;
    std::string read_1, read_2;
    int i = 1;
    gzFile fastain_1, fastain_2;

    fastain_1 = gzopen(rna_seq_1.c_str(), "r");
    seqs_1 = kseq_init(fastain_1);
    fastain_2 = gzopen(rna_seq_2.c_str(), "r");
    seqs_2 = kseq_init(fastain_2);

    int count1, count2;
    char strand1, strand2;

    // - Main loop: one iteration, one read
    // ---------------------------------------
    // l1 and l2 could probably be removed
    while ( ((l1 = kseq_read(seqs_1)) >= 0) && ((l2 = kseq_read(seqs_2)) >= 0) ) {

        count1 = 0;
        count2 = 0;

        head_1 = seqs_1->name.s;
        read_1 = seqs_1->seq.s;
        head_2 = seqs_2->name.s;
        read_2 = seqs_2->seq.s;

        strand1 = '+';
        strand2 = '+';

        // SEE: https://salmon.readthedocs.io/en/latest/library_type.html
        if (stranded) {
            if (inward && forward) {
                read_2 = reverseAndComplement(read_2);
                strand2 = '-';
            } else if (inward && reverse) {
                read_1 = reverseAndComplement(read_1);
                strand1 = '-';
            //} else if (matching && foward) {
                // do nothing
            } else if (matching && reverse) {
                read_1 = reverseAndComplement(read_1);
                strand1 = '-';
                read_2 = reverseAndComplement(read_2);
                strand2 = '-';
            } else if (outward && reverse) {
                read_1 = reverseAndComplement(read_1);
                strand1 = '-';
            } else if (outward && forward) {
                read_2 = reverseAndComplement(read_2);
                strand2 = '-';
            }

            // Only the read has to be aligned (as-is, see the previous if/else statements)
            paths = analyzeReadSimple(bm, sg, read_1, L, eps, exsN, verbose, strand1);

        } else if (unstranded) {
            // Both the read and its reverse-complement have to be aligned
            // in order to determine the best one
            paths = analyzeRead(bm, sg, read_1, L, eps, exsN, verbose);
        }

        // - Align first read
        // --------------------------------------------------------
        if(paths.first != '/') {
            for(std::pair<int, std::list<Mem> > path : paths.second) {
                count1++;
                if(!path.second.empty()) {
                    int err = path.first;
                    outFile_1 << "MAPPED" << " " << paths.first << " " << head_1 << " " << err << " ";
                    for(std::list<Mem>::iterator m=path.second.begin(); m!=path.second.end(); ++m) {
                        outFile_1 << m->toStr() << " ";
                    }

                    if (unstranded) {
                        if(paths.first == '+')
                            outFile_1 << read_1;
                        else
                            outFile_1 << reverseAndComplement(read_1);
                    } else {
                        // In case of a stranded library, the eventual reverse-and-complement
                        // has already been done, doesn't have to be done again
                        outFile_1 << read_1;
                    }

                    outFile_1 << "\n";
                }
            }
        } else {
            count1++;
            outFile_1 << "UNMAPPED" << " " << head_1 << " " << read_1 << "\n";
        }
        if(i%100 == 0)
            std::cout << "Processed " << i << " genes." << std::endl;
        ++i;

        if (unstranded) {
            if (inward) {
                if (paths.first == '+') {
                    read_2 = reverseAndComplement(read_2);
                    strand2 = '-';
                } //else if (paths.first == '-') {
                //} do nothing
            } else if (matching) {
                //if(paths.first == '+') {
                // do nothing
                //} else
                if (paths.first == '-') {
                    read_2 = reverseAndComplement(read_2);
                    strand2 = '-';
                }
            } else if (outward) {
                if (paths.first == '+') {
                    read_2 = reverseAndComplement(read_2);
                    strand2 = '-';
                } //else if (paths.first == '-') {
                  //do nothing
                  //}
            }
        }

        // - Align second read
        // --------------------------------------------------------

        // NOTE: in both stranded and unstranded cases the strand orientation of the second
        // read is already known
        paths = analyzeReadSimple(bm, sg, read_2, L, eps, exsN, verbose, strand2);

        if(paths.first != '/') {
            for(std::pair<int, std::list<Mem> > path : paths.second) {
                count2++;
                if(!path.second.empty()) {
                    int err = path.first;
                    outFile_2 << "MAPPED"<< " " << paths.first << " " << head_2 << " " << err << " ";
                    for(std::list<Mem>::iterator m=path.second.begin(); m!=path.second.end(); ++m) {
                        outFile_2 << m->toStr() << " ";
                    }

                    // Since the eventual reverse-and-complement has already been done in
                    // both stranded and unstranded cases, there is no need to do it again
                    outFile_2 << read_2 << "\n";

                } else {
                    int err = path.first;
                    outFile_2 << paths.first << " " << head_2 << " " << err << " ";
                    if(paths.first == '+')
                        outFile_2 << read_2;
                    else
                        outFile_2 << reverseAndComplement(read_2);
                }
            }
        } else {
            count2++;
            outFile_2 << "UNMAPPED" << " " << head_2 << " " << read_2 << "\n";
        }
        if(i%100 == 0)
            std::cout << "Processed " << i << " genes." << std::endl;
        ++i;

        if (count1 < count2) {
            while (count1 < count2) {
                outFile_1 << "PLACEHOLDER" << " " << head_1 << "\n";
                count1++;
            }
        } else if (count2 < count1) {
            while (count2 < count1) {
                    outFile_2 << "PLACEHOLDER" << " " << head_2 << "\n";
                    count2++;
            }
        }
    }

    kseq_destroy(seqs_1);
    gzclose(fastain_1);
    kseq_destroy(seqs_2);
    gzclose(fastain_2);

    return 0;
}
