# SGAL

Splicing Graph ALigner (SGAL) performs an accurate approximate mapping of RNA-seq reads against a graph representing the transcripts of a gene.

The version of the tool submitted to AlCoB 2017 can be found [here](https://github.com/AlgoLab/galig/tree/v1.0) (release v1.0).

### Prerequisites
  * python or python3
  * biopython

### Compiling
```bash
git clone https://github.com/AlgoLab/galig.git
cd galig
make prerequisites
make
```
### Running
```bash
# Pattern matching
./bin/main -g ./example/genomic.fa -a ./example/annotation.gtf -r ./example/rna_seqs.fa -l 3 -e 10 -o OUT
# SAM formatting
python3 scripts/SAMFormatter.py OUT ./example/genomic.fa.sg ./example/rna_seqs.fa
```

The tool has been tested only on 64bit Linux system.

[![Join the chat at https://gitter.im/AlgoLab/galig](https://badges.gitter.im/AlgoLab/galig.svg)](https://gitter.im/AlgoLab/galig?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
