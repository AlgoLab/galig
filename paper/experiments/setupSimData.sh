#!/bin/bash

WF=$1
DataFold=${WF}/SimData

log="${WF}/setupSimData.log"

# Reference
mkdir -p ${DataFold}/Refs

echo "*** Splitting reference"
python3 ./DownloadSimData/splitReference.py ${DataFold}/GRCh37.p13.genome.fa ${DataFold}/Refs 2> ${log}
rm ${DataFold}/GRCh37.p13.genome.fa

echo ""

# Annotation
mkdir -p ${DataFold}/Annotations

echo "*** Splitting annotation"
python3 ./DownloadSimData/splitAnnotation.py ${DataFold}/gencode.v19.annotation.mult_iso_subsample_1000_genes.gtf ${DataFold}/Annotations/ 2> ${log}
rm -f ${DataFold}/gencode.v19.annotation.mult_iso_subsample_1000_genes.gtf*

echo ""

# Samples
echo "*** Cleaning sample 5M"
bash ./DownloadSimData/cleanAndSplitSample.sh ${DataFold}/5000000_reads.fastq ${DataFold}/Samples/5M ${log}
rm -f ${DataFold}/5000000_reads.fasta &> ${log}

echo "*** Cleaning sample 10M"
bash ./DownloadSimData/cleanAndSplitSample.sh ${DataFold}/10000000_reads.fastq ${DataFold}/Samples/10M ${log}
rm -f ${DataFold}/10000000_reads.fasta &> ${log}

echo ""
