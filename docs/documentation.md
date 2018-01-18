[Main Page](index)

# Documentation

ASGAL has been developed and tested on Ubuntu Linux, but should work
on any *nix-based system.

<!--
* [Installation](#installation)
* [Usage](#usage)
* [Example](#example)
-->

## Installation
To compile ASGAL and the 3rd party libraries it uses
([sdsl-lite](https://github.com/simongog/sdsl-lite) and
[lemon](http://lemon.cs.elte.hu/trac/lemon)), we have to install:
  * [python3](https://www.python.org)
  * [gffutils](http://daler.github.io/gffutils/)
  * [biopython](http://biopython.org)
  * [cmake](https://cmake.org)

On an ubuntu system, something like this should be sufficient:
```bash
sudo apt-get update
sudo apt-get install build-essential cmake python3 python3-pip python3-biopython python3-biopython-sql
pip3 install --user gffutils
```

Then, to download and compile sdsl-lite, lemon and ASGAL:
```bash
git clone https://github.com/AlgoLab/galig.git
cd galig
make prerequisites
make
```
<br />

## Input
ASGAL takes as input:
* a reference genome (in FASTA format)
* a gene annotation (in GTF format)
* an RNA-Seq sample (in FASTA format)

##### FASTQ Support
If your sample is in FASTQ format, it can be easily converted into FASTA format using [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/):
```bash
fastq_to_fasta -i sample.fastq -o sample.fasta
```
<br />

## Usage
ASGAL is composed of different modules, each one performing a
different task.

### Step 1 - Splice-Aware Aligner

First of all, we have to build the splicing graph and align the input
sample to it:
```bash
./bin/main -g reference -a annotation -r sample -o output.mem
```

In this way, we obtain a file _output.mem_ containing the alignments
to the splicing graph.

### Step 2 (optional) - SAM Formatter
The alignments to the splicing graph are stored in a (pretty
unreadable) file. This file can be converted into a
[SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file (that can
be open, for example, with
[IGV](http://software.broadinstitute.org/software/igv/)) using:

```bash
python3 ./scripts/formatSAM.py output.mem annotation output.sam
```

##### Observation

The SAM file we produce contains the alignments to the splicing graph
mapped to the reference genome. These alignments must not be confused
with the spliced alignments to the reference genome computed by any
spliced aligner (such as [STAR](https://github.com/alexdobin/STAR) or
[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)).

Indeed, our alignments could contain long insertions representing the
portions of reads that could align to an intron (and represents a
cassette exon).

### Step 3 - Events Detector
Finally, we can analyze the alignments to detect the possible presence
of alternative splicing events:
```
python3 ./scripts/detectEvents.py reference annotation output.mem output.events
```

In this way, we obtain the file _output.events_ that contains the
detected alternative splicing events.

##### Output Format
Each alternative splicing event
is associated to the intron inducing it and it is described by:

 1. type of event:
   * **ES**, Exon Skipping
   * **A3**, Alternative acceptor (3') site
   * **A5**, Alternative donor (5') site
   * **IR**, Intron Retention
 2. genomic positions (start and end) of the intron inducing the event
 3. number of reads confirming the event
 4. list of annotated transcripts with respect to which the event is expressed

The _output.events_ file is a  _space-separated value_ file where each
line represents an event and each column represents one of the
previously describe features.

For example,

| Type | Start | End | #Reads | Annotated Transcripts |
|:-:|:-:|:-:|:-:|:-:|
| ES | 2000 | 2500 | 150 | Trancript1/Trancript2/Transcript3 |

this entry describes an exon skipping event, induced by 150 reads that
support the intron starting at position 2000 and ending at position 2500.
The event is expressed with respect to the annotated transcript
Trancript1, Trancript2, and Transcript3.

<br />

## Example
```bash
mkdir example
cd example
wget ftp://ftp.ensembl.org/pub/release-91/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.chromosome.X.fa.gz
gunzip Drosophila_melanogaster.BDGP6.dna.chromosome.X.fa.gz
mv Drosophila_melanogaster.BDGP6.dna.chromosome.X.fa DrosMel.chrX.fa
wget ftp://ftp.ensembl.org/pub/release-91/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.91.chr.gtf.gz
gunzip Drosophila_melanogaster.BDGP6.91.chr.gtf.gz
grep "CG13375" Drosophila_melanogaster.BDGP6.91.chr.gtf > CG13375.gtf
rm Drosophila_melanogaster.BDGP6.91.chr.gtf

../bin/main -g DrosMel.chrX.fa -a CG13375.gtf -s CG13375.fasta -o CG13375.mem
python3 ../scripts/formatSAM.py -m CG13375.mem -g DrosMel.chrX.fa -a CG13375.gtf -o CG13375.sam
python3 ../scripts/detectEvents.py -g DrosMel.chrX.fa -a CG13375.gtf -m CG13375.mem -o CG13375.events
```
<br />

## Command Line Arguments
* all programs show usage information with **-h** (**\-\-help**)
* all programs write to a specific output file, defined with **-o**

### Splice-Aware Aligner
Required parameters:
```
-g,--genome INFILE          FASTA input file containing the reference
-a,--annotation INFILE      GTF input file containing the gene annotation
-s,--sample INFILE          FASTA input file containing the RNA-Seq reads
-o,--output OUTFILE         output file containing the alignments to the
                            splicing graph
```

Optional parameters:
```
-l,--L <int>                minimum lenght of MEMs used to build the
                            alignments (default: 15)
-e,--eps <int>              error rate, a value from 0 to 100 used to
                            compute the maximum number of allowed errors
                            (default: 3)
```

##### Observations
* an higher value of _L_ improves the speed but reduces the number of aligned reads
* the number of allowed errors is computed for each read as the ratio between the error rate and its length

### SAM Formatter
Required parameters:
```
-g,--genome INFILE          FASTA input file containing the reference
-a,--annotation INFILE      GTF input file containing the gene annotation
-m,--mems INFILE            input file containing the alignments to
                            the splicing graph
-o,--output OUTFILE         SAM output file
```

Optional paramter:
```
-e,--erate <int>            error rate, a value from 0 to 100 (default: 3)
```
##### Observations
* the error rate should be the same used to compute the alignments

### Events Detector
Required parameters:
```
-g,--genome INFILE          FASTA input file containing the reference
-a,--annotation INFILE      GTF input file containing the gene annotation
-m,--mems INFILE            input file containing the alignments to
                            the splicing graph
-o,--output OUTFILE         output file containing the events
```

Optional parameters:
```
-e,--erate <int>            error rate, a value from 0 to 100 (default: 3)
-w,--support <int>          minimum number of reads needed to confirm an event
                            (default: 3)
--allevents                output all events, not only the novel ones
```

##### Observation
* the error rate should be the same used to compute the alignments
* by default, _asgal_ outputs only the alternative splicing events that
are **novel** with respect to the given annotation, that are the the
events induced by an intron not contained in the annotation