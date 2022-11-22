#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <utility>
#include <list>
#include <chrono>
#include <filesystem>

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
  std::string annotation;
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
	  {"genomic", required_argument, 0, 'g'},
	  {"annotation", required_argument, 0, 'a'},
	  {"annotation_dir", required_argument, 0, 'd'},
	  {"sample", required_argument, 0, 's'},
	  {"threads", required_argument, 0, 'j'},
	  {"L", required_argument, 0, 'l'},
	  {"erate", required_argument, 0, 'e'},
	  //{"output", required_argument, 0, 'o'},
	  {"help", no_argument, 0, 'h'},
	  //{"verbose", no_argument, 0, 'v'},
	  {0, 0, 0, 0}};

      int option_index = 0;
      c = getopt_long(argc, argv, "g:s:a:d:j:l:e:h", long_options, &option_index);

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
	  annotation = optarg;
	  break;
        case 'd':
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
  // std::string out_dir = annotation_dir+"/ASGAL/";

  // std::cout << "Genomic directory: " << genomic_dir << std::endl;
  // std::cout << "RNA-seq directory: " << rna_seqs_dir << std::endl;
  // std::cout << "Annotation file: " << annotation << std::endl;
  // std::cout << "Annotation directory: " << annotation_dir << std::endl;
  // // std::cout << "Output directory: " << out_dir << std::endl;

  // std::cout << "Minimum MEM length: " << L << std::endl;
  // std::cout << "Error rate: " << eps << std::endl;
  std::vector<SplicingGraph> sgs;
  std::vector<BackwardMEM> bms;
  std::map<std::string, std::string> chr_to_sg;
  // TODO parallelize by creating a class for sgs and bms
  // int count = 0;
  for (const auto &entry : std::filesystem::directory_iterator(genomic_dir))
    {
      // std::cout << entry.path() << std::endl;

      gzFile fastain = gzopen(entry.path().string().c_str(), "r");
      kseq_t *reference = kseq_init(fastain);
      kseq_read(reference);
      chr_to_sg[entry.path().stem().string()] = entry.path().string();
      kseq_destroy(reference);
      gzclose(fastain);
    }
  // std::cout << "Splicing graphs built" << std::endl;
  // std::cout << "Number of splicing graphs: " << sgs.size() << std::endl;
  // std::vector<std::pair<std::string, std::string>> rna_seqs;
  // for (const auto &entry : std::filesystem::directory_iterator(rna_seqs_dir))
  // {
  //     if (entry.path().extension().string() == ".fq")
  //     {
  //         rna_seqs.emplace_back(std::make_pair(entry.path().stem().string(), entry.path().string()));
  //     }
  // }
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
  if (verbose){
    for (auto &r : rna_seqs)
    {
      std::cout << r.name << "\t" << r.size << std::endl;
    }
    std::cout << "OMP thread max: " << omp_get_max_threads() << std::endl;
    std::cout << "thread max: " << std::thread::hardware_concurrency() << std::endl;
  }
  
  auto compute_mem = [annotation_dir, genomic_dir, &chr_to_sg, L, eps, verbose](rna_seq_read &sample)
  {
    if(verbose){
      std::cout << "sample: " << sample.path << std::endl;
    }
    auto start = std::chrono::high_resolution_clock::now();
    auto chr = chr_to_sg[split(sample.name, "_")[2]];
    std::string ann_tmp = annotation_dir + "/" + sample.name + "/annotation.gtf";
    std::string genome = split(sample.name, "_")[2];
    std::string genomic = genomic_dir + "/" + genome + ".fa";
   
    gzFile fastain = gzopen(genomic.c_str(), "r");
    kseq_t *reference = kseq_init(fastain);
    kseq_read(reference);

    SplicingGraph sg (reference->seq.s, ann_tmp);
    BackwardMEM bm(sg.getText(), genomic);
    //sg.print();
    kseq_destroy(reference);
    gzclose(fastain);

    int exsN = sg.getExonsNumber();

    gzFile fastain_q = gzopen(sample.path.c_str(), "r");
    //fastain = gzopen(sample.path.c_str(), "r");
    std::pair<char, std::list<std::pair<int, std::list<Mem>>>> paths;

    kseq_t *seqs = kseq_init(fastain_q);
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
        
    while ((l = kseq_read(seqs)) >= 0)
      {

	head = seqs->name.s;
	read = seqs->seq.s;
	paths = analyzeRead(bm, sg, read, L, eps, exsN, verbose);
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
    outFile.close();
    kseq_destroy(seqs);
    gzclose(fastain_q);
    
    std::string sam_command = "python3 scripts/formatSAM.py -m " + out + " -g "
      + chr + " -a " + ann_tmp + " -o "
      + annotation_dir + "/"
      + sample.name + "/ASGAL/aligns.sam";
    std::system(sam_command.c_str());
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    if (verbose){
      std::cout << "sample: " << sample.path << " time: " << duration.count() << std::endl;
    }
  };

  /* openmp */
  // #pragma omp parallel for
  //   for (auto &sample : rna_seqs)
  //     {
  //       compute_mem(std::ref(sample));
  //     }

  /* for_each c++17*/
  // // std::for_each(std::execution::par_unseq, rna_seqs.begin(), rna_seqs.end(), compute_mem);

  /* blocks of threads */
  // int n = (int) std::thread::hardware_concurrency()/2;
  // //int n = 200;
  // int size = (rna_seqs.size() - 1) / n + 1;
  // std::vector<rna_seq_read> vec[size];
  
  // // and store it in a vector at k'th index in `vec`
  // for (int k = 0; k < size; ++k)
  //   {
  //     // get range for the next set of `n` elements
  //     auto start_itr = std::next(rna_seqs.cbegin(), k*n);
  //     auto end_itr = std::next(rna_seqs.cbegin(), k*n + n);
      
  //     // allocate memory for the sub-vector
  //     vec[k].resize(n);
 
  //     // code to handle the last sub-vector as it might
  //     // contain fewer elements
  //     if (k*n + n > (int)rna_seqs.size())
  // 	{
  // 	  end_itr = rna_seqs.cend();
  // 	  vec[k].resize(rna_seqs.size() - k*n);
  // 	}
 
  //     // copy elements from the input range to the sub-vector
  //     std::copy(start_itr, end_itr, vec[k].begin());
  //   }
  // std::cout << rna_seqs.size() << " " << n << " " << size << "\n";
  // int ch = 0;
  // for (auto v: vec){
  //   std::cout << "pool number " << ch << "\n";
  //   std::vector<std::thread> th_vec;
  //   for (auto &sample : v)
  //     {
  // 	th_vec.emplace_back(compute_mem, std::ref(sample));
  //     }

  //   for (auto &th : th_vec)
  //     {
  // 	if(th.joinable())
  // 	  th.join();
  //     }
  //   ch++;
  //   th_vec.clear();
  // }

  //   std::vector<std::thread> th_vec;
  // for (auto &sample : rna_seqs)
  // {
  //     th_vec.emplace_back(compute_mem, std::ref(sample));
  // }

  // for (auto &th : th_vec)
  // {
  //     th.join();
  // }

 
  /* cptl */
  ctpl::thread_pool p(num_threads/2 - 1);
  for (auto &sample : rna_seqs)
    {
      p.push(std::bind(compute_mem, std::ref(sample)));
    }
  return 0;
}
