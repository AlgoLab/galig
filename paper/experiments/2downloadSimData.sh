#!/bin/bash

WF=$1
DataFold=${WF}/SimData

log="./getSimData.log"

# Reference
echo "*** Reference"
mkdir -p ${DataFold}/Refs

echo "** Downloading reference"
wget http://public.bmi.inf.ethz.ch/projects/2015/spladder/simulation/genome/GRCh37.p13.genome.fa.gz &> ${log}
echo "** Extracting reference"
gunzip -f GRCh37.p13.genome.fa.gz &> ${log}
echo "** Splitting reference"
python3 ./DownloadSimData/splitReference.py GRCh37.p13.genome.fa ${DataFold}/Refs 2> ${log}
rm GRCh37.p13.genome.fa

echo ""

# Annotation
echo "*** Annotation"
mkdir -p ${DataFold}/Annotations

echo "** Downloading annotation"
wget http://public.bmi.inf.ethz.ch/projects/2015/spladder/annotation/gencode.v19.annotation.mult_iso_subsample_1000_genes.gtf &> ${log}
echo "** Splitting annotation"
python3 ./DownloadSimData/splitAnnotation.py gencode.v19.annotation.mult_iso_subsample_1000_genes.gtf ${DataFold}/Annotations/ 2> ${log}
rm -f gencode.v19.annotation.mult_iso_subsample_1000_genes.gtf*

echo ""

# Samples
echo "*** Samples"
echo "** Downloading sample 5M"
wget http://public.bmi.inf.ethz.ch/projects/2015/spladder/simulation/reads/5000000_reads.fastq.gz &> ${log}
echo "** Extracting sample 5M"
gunzip -f 5000000_reads.fastq.gz &> ${log}
echo "** Cleaning sample 5M"
bash ./DownloadSimData/cleanAndSplitSample.sh 5000000_reads.fastq ${DataFold}/Samples/5M ${log}
rm -f 5000000_reads.fasta

echo "** Downloading sample 10M"
wget http://public.bmi.inf.ethz.ch/projects/2015/spladder/simulation/reads/10000000_reads.fastq.gz &> ${log}
echo "** Extracting sample 10M"
gunzip -f 10000000_reads.fastq.gz &> ${log}
echo "** Cleaning sample 10M"
bash ./DownloadSimData/cleanAndSplitSample.sh 10000000_reads.fastq ${DataFold}/Samples/10M ${log}
rm -f 10000000_reads.fasta &> ${log}

echo ""
