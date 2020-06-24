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
  * [python3](https://www.python.org) (>=3.6)
  * [biopython](http://biopython.org)
  * [pysam](https://pysam.readthedocs.io/en/latest/index.html)
  * [gffutils](http://daler.github.io/gffutils/)
  * [pandas](https://pandas.pydata.org/)
  * [cmake](https://cmake.org)
  * [samtools](http://samtools.sourceforge.net/)
  * [zlib](http://zlib.net/)

See <a href="http://asgal.algolab.eu/documentation#installation" target="_blank">here</a> for more details.

### Compiling
```bash
git clone --recursive https://github.com/AlgoLab/galig.git
cd galig
make prerequisites
make
```

### Running
```bash
./asgal -g [genome] -a [annotation] -s [sample] -o outputFolder
```

In more detail:
```bash
# Align RNA-Seq reads to a splicing graph
./bin/SpliceAwareAligner -g [reference] -a [annotation] -s [sample] -o outputFolder/output.mem

# Convert alignments to SAM format
python3 ./scripts/formatSAM.py -m output.mem -g [reference] -a [anotation] -o outputFolder/output.sam

# Detect events from alignments
python3 ./scripts/detectEvents.py -g [reference] -a [annotation] -m output.mem -o outputFolder/output.events.csv
```

### Example
```bash
cd example
tar xfz input.tar.gz
../asgal -g ./input/genome.fa -a ./input/annotation.gtf -s ./input/sample_1.fa -o outputFolder
```

This command will produce four files in the output folder:
  * a _.mem_ file containing the alignments to the splicing graph
  * a _.sam_ file, containing the alignments to the splicing graph mapped to the reference genome
  * a _.events.csv_ file, containing the alternative splicing events detected in the RNA-Seq sample
  * a _.log_ file, containing the log of the execution

An extended explanation of this example can be found <a href="http://asgal.algolab.eu/documentation#example" target="_blank">here</a>.


The tool has been tested only on 64bit Linux system. You can find more information at [http://asgal.algolab.eu](http://asgal.algolab.eu).

[![Join the chat at https://gitter.im/AlgoLab/galig](https://badges.gitter.im/AlgoLab/galig.svg)](https://gitter.im/AlgoLab/galig?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
