[//]: # (Comment)
# Genome-Wide Analysis

ASGAL tool is specifically designed to perform AS prediction based on
a splice-aware alignment of a RNA-seq sample against the splicing
graph of a specific gene. The current version of _ASGAL_ performs
efficiently in time when a limited set of genes are analyzed, while
for genome-wide analysis we have implemented a pre-processing step
that aims to speed up the process of filtering reads that map to genes
under investigation.

The genome-wide pipeline of ASGAL is mainly based on the **quasi-mapping**
algorithm of [Salmon](https://combine-lab.github.io/salmon/) and it
can be summarized as follows:

1. using the quasi-mapping, it quantifies the transcripts of the genes
and quick assigns each read to them

2. using the alignments produced by Salmon, it builds a set of
smaller samples, one for each considered gene

3. it splits the input references and annotations in multiple references and annotations

4. it runs ASGAL on each gene 

This pipeline can be run using the _ASGAL_GW_ script. Currently, this
script assumes that you have installed Salmon systemwide.

```bash
# Single-end sample
./ASGAL_GW -g [genome.fasta] \
           -a [annotation.gtf] \
           -s1 [sample.fa] \
           -t [transcripts.fasta] \
           -o [output_folder]

# Paired-end sample 
./ASGAL_GW -g [genome.fasta] \
           -a [annotation.gtf] \
           -s1 [sample1.fa] \
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

File:
```bash
./ASGAL_GW
```
Parameters:
```bash
-g,--genome INFILE          FASTA input file containing the chromosome sequences
-a,--annotation INFILE      GTF input file containing the annotations of considered genes
-s1,--sample1 INFILE        FASTA input file containing the RNA-Seq reads (sample 1)
-s2,--sample2 INFILE        FASTA input file containing the RNA-Seq reads (sample 2)
-t,--transcripts INFILE     FASTA input file containing the transcript sequences
                            of the considered genes
-o,--output OUTFOLD         output name folder
-l,--L <int>                minimum lenght of MEMs used to build the
                            alignments (default: 15)
-e,--eps <int>              error rate, a value from 0 to 100 used to
                            compute the maximum number of allowed errors
                            (default: 3)
-w,--support <int>          minimum number of reads needed to confirm an event
                            (default: 3)
-f,--allevents              output all events, not only the novel ones
                            (default: only novels)
```