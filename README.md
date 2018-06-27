# ASGAL

**ASGAL** (**A**lternative **S**plicing **G**raph **AL**igner) is a
tool for detecting the alternative splicing events expressed in a
RNA-Seq sample with respect to a gene annotation. The **main idea**
behind _ASGAL_ is the following one: the alternative splicing events can
be detected by aligning the RNA-Seq reads against the splicing graph
of the gene.

The instructions to install and use _ASGAL_ are at
[http://asgal.algolab.eu](http://asgal.algolab.eu).


### Prerequisites
  * [git](https://git-scm.com/)
  * [python3](https://www.python.org)
  * [gffutils](http://daler.github.io/gffutils/)
  * [biopython](http://biopython.org)
  * [cmake](https://cmake.org)
  * [zlib1g-dev](http://zlib.net/)

### Compiling
```bash
git clone --recursive https://github.com/AlgoLab/galig.git
cd galig
make prerequisites
make
```

### Running
```bash
./asgal -g [genome] -a [annotation] -s [sample] -o output
```

In more detail:
```bash
# Align RNA-Seq reads to a splicing graph
./bin/SpliceAwareAligner -g [reference] -a [annotation] -s [sample] -o output.mem

# Convert alignments to SAM format
python3 ./scripts/formatSAM.py -m output.mem -g [reference] -a [anotation] -o output.sam

# Detect events from alignments
python3 ./scripts/detectEvents.py -g [reference] -a [annotation] -m output.mem -o output.events.csv
```

### Example
```bash
cd example
tar xfz input.tar.gz
../asgal -g ./input/genome.fa -a ./input/annotation.gtf -s ./input/reads.fasta -o output
```

This command will produce three files:
  * _output.mem_, containing the alignments to the splicing graph
  * _output.sam_, containing the alignments to the splicing graph mapped to the reference genome
  * _output.events_, containing the alternative splicing events detected in the RNA-Seq sample

An extended explanation of this example can be found <a href="http://asgal.algolab.eu/documentation#example" target="_blank">here</a>.


The tool has been tested only on 64bit Linux system. You can find more information at [http://asgal.algolab.eu](http://asgal.algolab.eu).

[![Join the chat at https://gitter.im/AlgoLab/galig](https://badges.gitter.im/AlgoLab/galig.svg)](https://gitter.im/AlgoLab/galig?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
