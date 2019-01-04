[//]: # (Comment)
# Documentation
ASGAL is composed of different modules written in _c++_ and
_python_. It has been developed and tested on Ubuntu Linux, but should
work on any *nix-based system.

<!--
* [Installation](#installation)
* [Usage](#usage)
* [Example](#example)
-->

## Installation
To compile ASGAL and the 3rd party libraries it uses
([sdsl-lite](https://github.com/simongog/sdsl-lite) and
[lemon](http://lemon.cs.elte.hu/trac/lemon)), install:
  * [python3](https://www.python.org)
  * [biopython](http://biopython.org)
  * [pysam](https://pysam.readthedocs.io/en/latest/index.html)
  * [gffutils](http://daler.github.io/gffutils/)
  * [cmake](https://cmake.org)
  * [samtools](http://samtools.sourceforge.net/)
  * [zlib1g-dev](http://zlib.net/)

To compile [Salmon](https://combine-lab.github.io/salmon/) (used for
genome-wide analyses), install:
  * [boost](https://www.boost.org/)

On an ubuntu system (18.04), the following commands suffice:
```bash
sudo apt-get update
sudo apt-get install build-essential git python3 python3-pip python3-setuptools python3-biopython python3-biopython-sql python3-pysam cmake libboost1.65-all-dev samtools unzip wget curl zlib1g-dev liblzma-dev libjemalloc-dev libjemalloc1 libghc-bzlib-dev
pip3 install --user gffutils
```

Then, to download and compile sdsl-lite, lemon, Salmon and ASGAL:
```bash
git clone --recursive https://github.com/AlgoLab/galig.git
cd galig
make prerequisites
make
```

This creates the executable in the _bin_ folder. The python scripts,
instead, are in the _scripts_ folder.

<br />

## Input
ASGAL takes as input:
* a reference chromosome (in _FASTA_ format)
* a gene annotation (in _GTF_ format)
* an RNA-Seq sample (in _FASTA_ or _FASTQ_ format, it can be gzipped)

<br />

#### Note 1 (genome-wide analysis)
ASGAL tool takes as input the annotation of a single gene and the
relative chromosome: if you want to use the tool in a genome-wide
analysis, please refer to this [page](genomewide).

#### Note 2 (paired-end sample)
At the moment, ASGAL is not able to directly manage paired-end
samples: the only way to use both the fastq files (the reverse and the
forward one) consists in merging them into a unique fastq file and
then using this file as input for ASGAL. In this way, ASGAL will align
each read independently and then it will use the alignments to detect
the AS events. *Note that this is needed only if you run ASGAL on a
single gene: if you use the genome-wide pipeline, the merge is done
automatically.*

<br />

## Usage
ASGAL is composed of three different modules, each one performing a
different task. We made available a script that runs the full pipeline.

### Full Pipeline Script
To run ASGAL pipeline on a single gene, run the following command:
```bash
./asgal -g [genome] -a [annotation] -s [sample] -o outputFolder
```
This command will produce four files in the output folder:
  * a _.mem_ file containing the alignments to the splicing graph
  * a _.sam_ file, containing the alignments to the splicing graph mapped to the reference genome
  * a _.events.csv_ file, containing the alternative splicing events detected in the RNA-Seq sample
  * a _.log_ file, containing the log of the execution

We will now specify in more detail the three steps of ASGAL pipeline.

#### Step 1 - Splice-Aware Aligner
To build the splicing graph of the input gene and to align the input
sample to it, run the following command:
```bash
./bin/SpliceAwareAligner -g [reference] -a [annotation] -s [sample] -o outputFolder/output.mem
```

In this way, the alignments to the splicing graph are computed and
stored in the _output.mem_ file.

#### Step 2 (optional) - SAM Formatter
The file obtained in the previous step can be converted into a
[SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file (that can
be open, for example, with
[IGV](http://software.broadinstitute.org/software/igv/)) using:

```bash
python3 ./scripts/formatSAM.py -m output.mem -g [reference] -a [anotation] -o outputFolder/output.sam
```

###### Observation
The SAM file produce by ASGAL (_output.sam_) contains the alignments
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
python3 ./scripts/detectEvents.py -g [reference] -a [annotation] -m output.mem -o outputFolder/output.events.csv
```

###### Output Format
The alternative splicing events are stored in the file
_output.events.csv_, a _space-separated value_ file where each line
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

# Run ASGAL
../asgal -g input/genome.fa -a input/annotation.gtf -s input/reads.fasta -o CG13375
```

We should obtain a novel _exon skipping_ events:
```bash
Type,Start,End,Support,Transcripts
ES,287042,289040,411,FBtr0300326/FBtr0070103/FBtr0342963
```

###### Observation
The correctness of the output can be verified by opening the genome, the
annotation and the alignments with _IGV_ and by generating a _sashimi plot_:

<p align="center">
<img src="https://raw.githubusercontent.com/AlgoLab/galig/master/example/sashimi.png">
</p>

From the picture, we can see that 411 reads support the skipping of the third exon and this event involves the three transcript of the gene.

<br />

## Command Line Arguments
* all programs show usage information with **-h** (**\-\-help**)

### Full Pipeline Script

File:
```bash
./asgal
```
Required parameters:
```bash
-g,--genome INFILE          FASTA input file containing the reference
-a,--annotation INFILE      GTF input file containing the gene annotation
-s,--sample INFILE          FASTA/Q input file containing the RNA-Seq reads (can be a gzipped file)
-o,--output OUTFOLD         output folder
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
-s,--sample INFILE          FASTA/Q input file containing the RNA-Seq reads (can be a gzipped file)
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
* an higher value of _L_ improves the speed but reduces the number of aligned reads
* the number of allowed errors is computed for each read as the ratio between the error rate and its length

<br />

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

<br />

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
                            (default: only novels)
```

###### Observation
* the error rate should be the same used to compute the alignments
* by default, ASGAL outputs only the alternative splicing events that
are **novel** with respect to the given annotation, that are the the
events induced by an intron not contained in the annotation
