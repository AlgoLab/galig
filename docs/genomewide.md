[//]: # (Comment)
# Genome-Wide Analysis

ASGAL tool is specifically designed to perform AS prediction based on
a splice-aware alignment of a RNA-seq sample against the splicing
graph of a specific gene. The current version of _ASGAL_ performs
efficiently in time when a limited set of genes are analyzed, while
for genome-wide analysis we have implemented a pre-processing step
that aims to speed up the process of filtering reads that map to genes
under investigation.

The genome-wide mode of ASGAL is mainly based on the **quasi-mapping**
algorithm of [Salmon](https://combine-lab.github.io/salmon/) and it
can be summarized as follows:

1. using the quasi-mapping, it quantifies the transcripts of the genes
and quick assigns each read to them

2. using the alignments produced by Salmon, it builds a set of
smaller samples, one for each considered gene

3. it splits the input references and annotations in multiple references and annotations

4. it runs ASGAL on each gene

ASGAL can be run in genome-wide mode by passing to the _asgal_ script (the same used to
run ASGAL on a single gene) the "--multi" flag:
```bash
# Single-end sample
./asgal --multi \
        -g [genome.fasta] \
        -a [annotation.gtf] \
        -s [sample.fa] \
        -t [transcripts.fasta] \
        -o [output_folder]

# Paired-end sample
./asgal --multi \
        -g [genome.fasta] \
        -a [annotation.gtf] \
        -s [sample1.fa] \
        -s2 [sample2.fa] \
        -t [transcripts.fasta] \
        -o [output_folder]
```

This script takes as input:
* the chromosome sequences
* the annotation of the considered genes
* the RNA-Seq samples
* the transcripts of the considered genes
* the output folder

and it produces in the output folder:
* the _refs_ folder which contains the chromosome sequences
* the _annos_ folder which contains the annotations, one for each considered gene
* the _samples_ folder which contains the samples, one for each considered gene
* the _salmon_ folder which contains the outputs of Salmon
* the _ASGAL_ folder which contains the outputs of ASGAL
* the _logs_ folder containing the logs file for the different steps of the pipeline

## Command Line Arguments
* the script shows usage information with **-h** (**\-\-help**)

#### Warning
The transcript IDs contained in the input annotation and in the input
transcripts fasta file should match. If not, the script won't split
the input sample in smaller samples and ASGAL won't produce any
output.

File:
```bash
./asgal
```
Parameters:
```bash
--multi                     use this to run ASGAL in genome-wide mode
-g,--genome INFILE          FASTA input file containing the chromosome sequences
-a,--annotation INFILE      GTF input file containing the annotations of considered genes
-s,--sample INFILE          FASTA input file containing the first RNA-Seq sample (can be gzipped)
-s2,--sample2 INFILE        FASTA input file containing the second RNA-Seq sample (can be gzipped)
-t,--transcripts INFILE     FASTA input file containing the transcript sequences of the
                            considered genes (can be gzipped)
-o,--output OUTFOLD         output name folder
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

## Example
We built a simple example using 19 genes from human chromosomes 13 and Y. To run the example:

1 download the example data from [here](https://drive.google.com/open?id=1CX_9-a0-vUa2pImRePHQhFqEywO9DzTu)

2 unzip the archive:
```bash
tar xfz GW_ASGAL_example.tar.gz
cd GW_ASGAL_example
```
3 run the ASGAL pipeline:
```bash
/path/to/asgal --multi -g genome.fa -a annotations.gtf -s sample1.fa.gz -t transcripts.fa.gz -o outFold
```