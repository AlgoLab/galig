# ASGAL

**ASGAL** (**A**lternative **S**plicing **G**raph **AL**igner) is a
tool for detecting the alternative splicing events expressed in a
RNA-Seq sample with respect to a gene annotation. The **main idea**
behind _ASGAL_ is the following one: the alternative splicing events can
be detected by aligning the RNA-Seq reads against the splicing graph
of the gene.

### Prerequisites
  * [git](https://git-scm.com/)
  * [python3](https://www.python.org)
  * [gffutils](http://daler.github.io/gffutils/)
  * [biopython](http://biopython.org)
  * [cmake](https://cmake.org)

### Compiling
```bash
git clone https://github.com/AlgoLab/galig.git
cd galig
make prerequisites
make
```

### Running
```bash
# Align RNA-Seq reads to a splicing graph
./bin/main -g [reference] -a [annotation] -r [sample] -o output.mem

# Convert alignments to SAM format
python3 ./scripts/formatSAM.py -m output.mem -g [reference] -a [anotation] -o output.sam

# Detect events from alignments
python3 ./scripts/detectEvents.py -g [reference] -a [annotation] -m output.mem -o output.events 
```

### Example
```bash
cd example
tar xfJ DrosMel.BDGP6.chrX.tar.xz
../bin/main -g DrosMel.BDGP6.chrX.fa -a CG13375.gtf -s CG13375.fasta -o CG13375.mem
python3 ../scripts/formatSAM.py -m CG13375.mem -g DrosMel.BDGP6.chrX.fa -a CG13375.gtf -o CG13375.sam
python3 ../scripts/detectEvents.py --allevents -g DrosMel.BDGP6.chrX.fa -a CG13375.gtf -m CG13375.mem -o CG13375.events
```

An extended explanation of this example can be found <a href="http://asgal.algolab.eu/documentation#example" target="_blank">here</a>.


The tool has been tested only on 64bit Linux system.

You can find more information at [http://asgal.algolab.eu](http://asgal.algolab.eu).

[![Join the chat at https://gitter.im/AlgoLab/galig](https://badges.gitter.im/AlgoLab/galig.svg)](https://gitter.im/AlgoLab/galig?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
