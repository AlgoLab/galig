# SGAL

Splicing Graph ALigner (SGAL) performs an accurate approximate mapping of RNA-seq reads against a graph representing the transcripts of a gene.

This version of the tool is the one submitted to AlCoB 2017 (http://grammars.grlmc.com/AlCoB2017/).

### Prerequisites
  * python3 (maybe also python is good enough but I didn't test it)
  * python3 modules: biopython, bcbio-gff, editdistance

### Compiling
```bash
git clone --recursive https://github.com/AlgoLab/galig.git
cd galig
make prerequisites
make
cd example/
./run genomic.fa annotation.gff rna_seqs.fa 3 5 gene out
```

The tool has been tested only on 64bit Linux system.

[![Join the chat at https://gitter.im/AlgoLab/galig](https://badges.gitter.im/AlgoLab/galig.svg)](https://gitter.im/AlgoLab/galig?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
