#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <utility>
#include <list>
#include <chrono>
#include <filesystem>
#include <mutex>
#include <zlib.h>
#include <stdio.h>
#include <functional>
#include "kseq.h"
#include "../ctpl/ctpl_stl.h"
#include "SplicingGraph.hpp"
#include "bMEM.hpp"
#include "utils.hpp"
#include "MEMsGraph.hpp"
#include "omp.h"


KSEQ_INIT(gzFile, gzread)

//std::mutex mu;

struct rna_seq_read
{
  std::string name;
  std::string path;
  double size;

  rna_seq_read(std::string name, std::string path, double size) : name(name), path(path), size(size) {}
  rna_seq_read()=default;
};

bool size_compare(const rna_seq_read &x, const rna_seq_read &y) { return x.size > y.size; }

std::vector<std::string> split(std::string s, std::string delimiter)
{
  size_t pos_start = 0, pos_end, delim_len = delimiter.length();
  std::string token;
  std::vector<std::string> res;

  while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos)
    {
      token = s.substr(pos_start, pos_end - pos_start);
      pos_start = pos_end + delim_len;
      res.push_back(token);
    }

  res.push_back(s.substr(pos_start));
  return res;
}

void printHelp()
{
  std::cout << "Usage: SGAL [options] (required: -g -a -s -o)\n"
	    << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -g, --genomes <path>" << std::endl;
  std::cout << "  -a, --annotation <path>" << std::endl;
  std::cout << "  -s, --samples <path>" << std::endl;
  // std::cout << "  -o, --output <path>: output file" << std::endl;
  std::cout << "  -j, --threads <int>: maximum threads per pool"  << std::endl;
  std::cout << "  -l, --L <int>: minimum lenght of MEMs used to build the alignments (default: 15)" << std::endl;
  std::cout << "  -e, --eps <int>: error rate, a value from 0 to 100 (default: 3)" << std::endl;
  std::cout << "  -h, --help: show this help message and exit" << std::endl;
  // std::cout << "  -v, --verbose: explain what is being done and save .dot" << std::endl;
}

std::pair<char, std::list<std::pair<int, std::list<Mem>>>> analyzeRead(BackwardMEM &bm,
                                                                       const SplicingGraph &sg,
                                                                       const std::string &read,
                                                                       const int &L,
                                                                       const int &eps,
                                                                       const int &exsN,
                                                                       const bool &verbose)
{
  // Original read
  std::list<Mem> mems = bm.getMEMs(read, L);
  if (verbose)
    {
      for (const Mem &m : mems)
        {
	  std::cout << m.toStr() << " ";
        }
      std::cout << std::endl;
    }
  std::list<std::pair<int, std::list<Mem>>> paths; // Path: [(weight, [mems])]
  if (!mems.empty())
    {
      MemsGraph mg(read, L, eps, exsN, verbose);
      mg.build(sg, mems);
      paths = mg.visit(sg);
    }
  // Reversed-and-complemented read
  std::string readRC = reverseAndComplement(read);
  std::list<Mem> memsRC = bm.getMEMs(readRC, L);
  if (verbose)
    {
      for (const Mem &m : memsRC)
        {
	  std::cout << m.toStr() << " ";
        }
      std::cout << std::endl;
    }
  std::list<std::pair<int, std::list<Mem>>> pathsRC; // Path: [(weight, [mems])]
  if (!memsRC.empty())
    {
      MemsGraph mgRC(readRC, L, eps, exsN, verbose);
      mgRC.build(sg, memsRC);
      pathsRC = mgRC.visit(sg);
    }
  bool empty = paths.empty();
  bool emptyRC = pathsRC.empty();
  char strand = '/';
  if (!empty || !emptyRC)
    {
      if (!empty && emptyRC)
        {
	  strand = '+';
        }
      else if (empty && !emptyRC)
        {
	  paths = pathsRC;
	  strand = '-';
        }
      else
        {
	  if (paths.front().first <= pathsRC.front().first)
            {
	      strand = '+';
            }
	  else
            {
	      paths = pathsRC;
	      strand = '-';
            }
        }
    }
  return std::make_pair(strand, paths);
}

int main(int argc, char *argv[])
{
  std::string genomic_dir;
  std::string annotation_dir;
  std::string rna_seqs_dir;
  int num_threads = 0;
  int L = 0;
  int eps = -1;

  bool verbose = false;

  // - Collecting command line parameters
  // -------------------------------------
  int c;
  while (1)
    {
      static struct option long_options[] =
	{
	  {"genomic_dir", required_argument, 0, 'g'},
	  {"annotation_dir", required_argument, 0, 'a'},
	  {"sample_dir", required_argument, 0, 's'},
	  {"threads", required_argument, 0, 'j'},
	  {"L", required_argument, 0, 'l'},
	  {"erate", required_argument, 0, 'e'},
	  {"help", no_argument, 0, 'h'},
	  //{"verbose", no_argument, 0, 'v'},
	  {0, 0, 0, 0}};

      int option_index = 0;
      c = getopt_long(argc, argv, "g:s:a:j:l:e:h", long_options, &option_index);

      if (c == -1)
        {
	  break;
        }

      switch (c)
        {
        case 'g':
	  genomic_dir = optarg;
	  break;
        case 's':
	  rna_seqs_dir = optarg;
	  break;
        case 'a':
	  annotation_dir = optarg;
	  break;
	case 'j':
	  num_threads = std::stoi(optarg);
	  break;
        case 'l':
	  L = std::stoi(optarg);
	  break;
        case 'e':
	  eps = std::stoi(optarg);
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
  if (L == 0)
    {
      L = 15;
    }
  if (eps < 0 || eps > 100)
    {
      eps = 3;
    }
  if (num_threads == 0)
    {
      num_threads = (std::thread::hardware_concurrency() / 2) - 1;
    }

  std::vector<rna_seq_read> rna_seqs;
  for (const auto &entry : std::filesystem::directory_iterator(rna_seqs_dir))
    {
      if (entry.path().extension().string() == ".fq")
        {
	  rna_seqs.emplace_back(rna_seq_read(entry.path().stem().string(),
					     entry.path().string(),
					     std::filesystem::file_size(entry.path().string())));
        }
    }
  std::sort(rna_seqs.begin(), rna_seqs.end(), size_compare);
  auto compute_mem = [annotation_dir, genomic_dir, L, eps, verbose](int id, rna_seq_read sample)
  {
   auto start = std::chrono::high_resolution_clock::now();
    std::string ann_tmp = annotation_dir + "/" + sample.name + "/annotation.gtf";
    std::string genome = split(sample.name, "_")[2];
    std::string genomic = "";
    if(genomic_dir.back() != '/'){
      genomic = genomic_dir + "/" + genome + ".fa";
    }else{
      genomic = genomic_dir + genome + ".fa";
    }
   
    gzFile fastain = gzopen(genomic.c_str(), "r");
    kseq_t *reference = kseq_init(fastain);
    kseq_read(reference);
    //mu.lock();
    SplicingGraph sg (reference->seq.s, ann_tmp);
    //mu.unlock();
    kseq_destroy(reference);
    gzclose(fastain);
    BackwardMEM bm(sg.getText(), genomic);

    gzFile fastain_q = gzopen(sample.path.c_str(), "r");
    std::pair<char, std::list<std::pair<int, std::list<Mem>>>> paths;
    int l;
    std::string head;
    std::string read;
    int i = 1;
    
    if (!std::filesystem::is_directory(annotation_dir + "/" + sample.name + "/ASGAL") ||
	!std::filesystem::exists(annotation_dir + "/" + sample.name + "/ASGAL"))
      {
	std::filesystem::create_directory(annotation_dir + "/" + sample.name + "/ASGAL");
      }
    
    std::string out = annotation_dir + "/" + sample.name + "/ASGAL/aligns.mem";

    std::ofstream outFile;
    outFile.open(out);
    std::vector<std::pair<std::string, std::string>> reads;    
    kseq_t *seqs = kseq_init(fastain_q); 
    while ((l = kseq_read(seqs)) >= 0)
      {

	head = seqs->name.s;
	read = seqs->seq.s;
	paths = analyzeRead(bm, sg, read, L, eps, 0, verbose);
	if (paths.first != '/')
	  {
	    for (std::pair<int, std::list<Mem>> path : paths.second)
	      {
		if (!path.second.empty())
		  {
		    int err = path.first;
		    outFile << paths.first << " " << head << " " << err << " ";
		    for (std::list<Mem>::iterator m = path.second.begin(); m != path.second.end(); ++m)
		      {
			outFile << m->toStr() << " ";
		      }
		    if (paths.first == '+')
		      outFile << read;
		    else
		      outFile << reverseAndComplement(read);
		    outFile << "\n";
		  }
	      }
	  }
	if (i % 10000 == 0)
	  std::cout << "Processed " << i << " reads." << std::endl;
	++i;
      }
    kseq_destroy(seqs);
    gzclose(fastain_q);
    
    outFile.close();
    std::string sam_command = "python3 scripts/formatSAM.py -m " + out + " -g "
      + genomic + " -a " + ann_tmp + " -o "
      + annotation_dir + "/"
      + sample.name + "/ASGAL/aligns.sam";
    //int sam_ret = std::system(sam_command.c_str());
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
     if (false){
      std::cout << "sample: " << sample.path << "\n genome: " << genomic <<
	"\n annotation:" << ann_tmp  << "\n time: " << duration.count() <<
	std::endl << std::endl;
    }
  };
  /* openmp */
  // #pragma omp parallel for
  //   for (auto &sample : rna_seqs)
  //     {
  //       compute_mem(std::ref(sample));
  //     }

  /* cptl */
  ctpl::thread_pool p(num_threads);
  for (auto &sample : rna_seqs)
    {
      p.push(compute_mem, sample);
    }
  return 0;
}
