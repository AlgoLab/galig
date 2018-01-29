[Main Page](index)

# Documentation
_ASGAL_ is composed of different modules written in _c++_ and
_python_. It has been developed and tested on Ubuntu Linux, but should
work on any *nix-based system.

<!--
* [Installation](#installation)
* [Usage](#usage)
* [Example](#example)
-->

## Installation
To compile _ASGAL_ and the 3rd party libraries it uses
([sdsl-lite](https://github.com/simongog/sdsl-lite) and
[lemon](http://lemon.cs.elte.hu/trac/lemon)), install:
  * [python3](https://www.python.org)
  * [gffutils](http://daler.github.io/gffutils/)
  * [biopython](http://biopython.org)
  * [cmake](https://cmake.org)

On an ubuntu system, the following commands suffice:
```bash
sudo apt-get update
sudo apt-get install build-essential cmake python3 python3-pip python3-biopython python3-biopython-sql
pip3 install --user gffutils
```

Then, to download and compile _sdsl-lite_, _lemon_ and _ASGAL_:
```bash
git clone https://github.com/AlgoLab/galig.git
cd galig
make prerequisites
make
```

This creates the executable in the _bin_ folder. The python scripts,
instead, are in the _scripts_ folder.

<br />

## Input
_ASGAL_ takes as input:
* a reference genome (in _FASTA_ or _FASTQ_ format)
* a gene annotation (in _GTF_ format)
* an RNA-Seq sample (in _FASTA_ format)

<br />

## Usage
_ASGAL_ is composed of three different modules, each one performing a
different task. However, we make available a script that run the full pipeline.

### Full Pipeline Script
To run _ASGAL_ pipeline, run the following command:
```bash
./asgal -g [genome] -a [annotation] -s [sample] -o output
```
This command will produce three files:
  * _output.mem_, containing the alignments to the splicing graph
  * _output.sam_, containing the alignments to the splicing graph mapped to the reference genome
  * _output.events.csv_, containing the alternative splicing events detected in the RNA-Seq sample

We will now specify in more detail the three steps of _ASGAL_ pipeline.

#### Step 1 - Splice-Aware Aligner

To build the splicing graph and align the input
sample to it, run the following command:
```bash
./bin/SpliceAwareAligner -g [reference] -a [annotation] -r [sample] -o output.mem
```

In this way, the alignments to the splicing graph are computed and
stored in _output.mem_ file.

###### Observation
The SpliceAwareAligner takes as input only FASTA samples. If you have
a FASTQ sample, you can convert it into FASTA format using the
following command:
```bash
python3 ./scripts/fastq2fasta.py -i [FASTQ] -o [FASTA]
```

#### Step 2 (optional) - SAM Formatter
The file obtained in the previous step can be converted into a
[SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file (that can
be open, for example, with
[IGV](http://software.broadinstitute.org/software/igv/)) using:

```bash
python3 ./scripts/formatSAM.py -m output.mem -g [reference] -a [anotation] -o output.sam
```

###### Observation
The SAM file produce by _ASGAL_ (_output.sam_) contains the alignments
to the splicing graph mapped to the reference genome. These alignments
must not be confused with the spliced alignments to the reference
genome computed by any spliced aligner (such as
[STAR](https://github.com/alexdobin/STAR) or
[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)).
Indeed, our alignments could contain long insertions representing the
portions of reads that could align to an intron (and that could
represent, for example, a _cassette exon_).

#### Step 3 - Events Detector

To analyze the alignments and detect the possible presence of
alternative splicing events, run:

```bash
python3 ./scripts/detectEvents.py -g [reference] -a [annotation] -m output.mem -o output.events.csv
```

###### Output Format
The alternative splicing events are stored in the file
_output.events_, a _space-separated value_ file where each line
represents an event. Each alternative splicing event is associated to
the intron inducing it and it is described by:
 1. type of event:
   * **ES**, Exon Skipping
   * **A3**, Alternative acceptor (3') site
   * **A5**, Alternative donor (5') site
   * **IR**, Intron Retention
 2. genomic positions (start and end) of the intron inducing the event
 3. number of reads confirming the event
 4. list of annotated transcripts with respect to which the event is expressed

For example,

| Type | Start | End | Support | Transcripts |
|:-:|:-:|:-:|:-:|:-:|
| ES | 2000 | 2500 | 150 | Transcript1/Transcript2/Transcript3 |

this entry describes an **exon skipping** event, induced by **150**
reads that support the intron starting at position **2000** and ending
at position **2500**.  The event is expressed with respect to the
annotated transcript **Transcript1, Transcript2, and Transcript3**.

<br />

## Example

We want align the RNA-Seq sample contained in the _example_ folder to
gene
[CG13375](http://www.ensembl.org/Drosophila_melanogaster/Gene/Summary?db=core;g=FBgn0040370;r=X:283186-294962)
of _Drosophila Melanogaster_.

```bash
cd example

# Extract the RNA-Seq sample
tar xfz input.tar.gz

# Download reference and annotation from ensembl
wget ftp://ftp.ensembl.org/pub/release-91/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.chromosome.X.fa.gz
gunzip Drosophila_melanogaster.BDGP6.dna.chromosome.X.fa.gz
mv Drosophila_melanogaster.BDGP6.dna.chromosome.X.fa input/genome.fa
wget ftp://ftp.ensembl.org/pub/release-91/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.91.chr.gtf.gz
gunzip Drosophila_melanogaster.BDGP6.91.chr.gtf.gz

# Extract the annotation of the gene of interest
grep "CG13375" Drosophila_melanogaster.BDGP6.91.chr.gtf > input/annotation.gtf
rm Drosophila_melanogaster.BDGP6.91.chr.gtf

# Run _ASGAL_
../bin/SpliceAwareAligner -g input/genome.fa -a input/annotation.gtf -s input/reads.fasta -o CG13375.mem
python3 ../scripts/formatSAM.py -m CG13375.mem -g input/genome.fa -a input/annotation.gtf -o CG13375.sam
python3 ../scripts/detectEvents.py -g input/genome.fa -a input/annotation.gtf -m CG13375.mem -o CG13375.events
```

We should obtain a novel _exon skipping_ events:
```bash
Type,Start,End,Support,Transcripts
ES,287042,289040,411,FBtr0300326/FBtr0070103/FBtr0342963
```

###### Observation
The correctness of the output can be verified opening the genome, the
annotation and the alignments with _IGV_ and generating a _sashimi plot_:

<img src="https://raw.githubusercontent.com/AlgoLab/galig/master/example/sashimi.png" width="630">

From the picture, we can see that 411 reads support the skipping of the third exon and this event involves the three transcript of the gene.

<br />

## Command Line Arguments
* all programs show usage information with **-h** (**\-\-help**)
* all programs write to a specific output file, defined with **-o**

### Full Pipeline script

File:
```bash
./asgal
```
Required parameters:
```bash
-g,--genome INFILE          FASTA input file containing the reference
-a,--annotation INFILE      GTF input file containing the gene annotation
-s,--sample INFILE          FASTA/Q input file containing the RNA-Seq reads
-o,--output OUTFILE         output file name (without extension)
```

Optional parameters:
```bash
-l,--L <int>                minimum lenght of MEMs used to build the
                            alignments (default: 15)
-e,--eps <int>              error rate, a value from 0 to 100 used to
                            compute the maximum number of allowed errors
                            (default: 3)
-w,--support <int>          minimum number of reads needed to confirm an event
                            (default: 3)
--allevents                 output all events, not only the novel ones
                            (default: only novels)
```

### Splice-Aware Aligner

File:
```bash
./bin/SpliceAwareAligner
```

Required parameters:
```bash
-g,--genome INFILE          FASTA input file containing the reference
-a,--annotation INFILE      GTF input file containing the gene annotation
-s,--sample INFILE          FASTA input file containing the RNA-Seq reads
-o,--output OUTFILE         output file containing the alignments to the
                            splicing graph
```

Optional parameters:
```bash
-l,--L <int>                minimum lenght of MEMs used to build the
                            alignments (default: 15)
-e,--eps <int>              error rate, a value from 0 to 100 used to
                            compute the maximum number of allowed errors
                            (default: 3)
```

###### Observations
* the SpliceAwareAligner takes as input only FASTA samples. If you have a FASTQ sample, you can convert it into FASTA format using the following command:
  ```bash
  python3 ./scripts/fastq2fasta.py -i [FASTQ] -o [FASTA]
  ```
* an higher value of _L_ improves the speed but reduces the number of aligned reads
* the number of allowed errors is computed for each read as the ratio between the error rate and its length

### SAM Formatter
File:
```bash
./scripts/formatSAM.py
```

Required parameters:
```bash
-g,--genome INFILE          FASTA input file containing the reference
-a,--annotation INFILE      GTF input file containing the gene annotation
-m,--mems INFILE            input file containing the alignments to
                            the splicing graph
-o,--output OUTFILE         SAM output file
```

Optional paramter:
```bash
-e,--erate <int>            error rate, a value from 0 to 100 (default: 3)
```
###### Observations
* the error rate should be the same used to compute the alignments

### Events Detector
File:
```bash
./scripts/detectEvents.py
```

Required parameters:
```bash
-g,--genome INFILE          FASTA input file containing the reference
-a,--annotation INFILE      GTF input file containing the gene annotation
-m,--mems INFILE            input file containing the alignments to
                            the splicing graph
-o,--output OUTFILE         output file containing the events
```

Optional parameters:
```bash
-e,--erate <int>            error rate, a value from 0 to 100 (default: 3)
-w,--support <int>          minimum number of reads needed to confirm an event
                            (default: 3)
--allevents                 output all events, not only the novel ones
```

###### Observation
* the error rate should be the same used to compute the alignments
* by default, _ASGAL_ outputs only the alternative splicing events that
are **novel** with respect to the given annotation, that are the the
events induced by an intron not contained in the annotation