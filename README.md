# SGAL

A Splicing Graph ALigner

Clone and compile with
```bash
git clone --recursive https://github.com/AlgoLab/galig.git
cd galig
make prerequisites
make
cd example/
./run genomic.fa annotation.gff rna_seqs.fa 3 5 gene out
```
or
```bash
git clone https://github.com/AlgoLab/galig.git
cd galig
git clone --recursive https://github.com/ldenti/backwardMEM.git
make prerequisites
make
cd example/
./run genomic.fa annotation.gff rna_seqs.fa 3 5 gene out
```

[![Join the chat at https://gitter.im/AlgoLab/galig](https://badges.gitter.im/AlgoLab/galig.svg)](https://gitter.im/AlgoLab/galig?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
