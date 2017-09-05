# SGAL

Splicing Graph ALigner (SGAL) performs an accurate approximate mapping of RNA-seq reads against a graph representing the transcripts of a gene. Moreover, it determines which alternative splicing events are expressed in the reads sample with respect to the known annotation of the gene. Finally, if more samples are given in input, the tool is able to compare the events extracted from each of them and outputs a summary table.

The version of the tool submitted to AlCoB 2017 can be found [here](https://github.com/AlgoLab/galig/tree/v1.0) (release v1.0).

### Prerequisites
  * python or python3
  * biopython
  * cmake

### Compiling
```bash
git clone https://github.com/AlgoLab/galig.git
cd galig
make prerequisites
make
```
### Running
```bash
./run -g ./example/genomic.fa -a ./example/annotation.gtf -r ./example/rna_seqs.fa -l 3 -e 10 -o OUT
```

More in detail...
```bash
# Pattern matching
./bin/main -g ./example/genomic.fa -a ./example/annotation.gtf -r ./example/rna_seqs.fa -l 3 -e 10 -o ./out
# SAM formatting
python3 ./scripts/SAMFormatter.py ./out ./example/annotation.gtf.sg ./example/rna_seqs.fa
# Events extractor
python3 ./scripts/postprocessing.py ./example/annotation.gtf.sg ./out > ./events
# Events comparator
python3 ./scripts/compareEvents.py ./events > comparison.csv
```

The tool has been tested only on 64bit Linux system.

[![Join the chat at https://gitter.im/AlgoLab/galig](https://badges.gitter.im/AlgoLab/galig.svg)](https://gitter.im/AlgoLab/galig?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)