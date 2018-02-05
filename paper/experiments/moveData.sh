#!/bin/bash

CF=$(pwd)

InFold=$1
OutFold=$2

SimFold=${OutFold}/SimData
RealFold=${OutFold}/RealData

mkdir -p ${SimFold}
mkdir -p ${RealFold}

cd ${SimFold}
mv ${InFold}/GRCh37.p13.genome.fa.gz .
gunzip -f GRCh37.p13.genome.fa.gz

mv ${InFold}/gencode.v19.annotation.mult_iso_subsample_1000_genes.gtf .

mv ${InFold}/5000000_reads.fastq.gz .
gunzip -f 5000000_reads.fastq.gz
mv ${InFold}/10000000_reads.fastq.gz .
gunzip -f 10000000_reads.fastq.gz

cd ${RealFold}

mv ${InFold}/SRR354041.tar.gz .
tar xfz SRR354041.tar.gz
rm SRR354041.tar.gz
mv ${InFold}/SRR354042.tar.gz .
tar xfz SRR354042.tar.gz
rm SRR354042.tar.gz
mv ${InFold}/SRR354043.tar.gz .
tar xfz SRR354043.tar.gz
rm SRR354043.tar.gz

cd ${CF}
