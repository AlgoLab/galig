#!/bin/bash

CF=$(pwd)

WF=$1

mkdir -p $WF
cd $WF

tar xfz genes_information.tar.gz

# Reference
mkdir -p Reference/Genomes
cd Reference/genomes
for chr in $(seq 1 22) X Y
do
    wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.${chr}.fa.gz
    gunzip Homo_sapiens.GRCh37.75.dna.chromosome.${chr}.fa.gz
done
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
cd ..

# Transcripts
wget https://github.com/comprna/SUPPA_supplementary_data/raw/master/annotation/hg19_EnsenmblGenes_sequence_ensenmbl.fasta.gz
gunzip hg19_EnsenmblGenes_sequence_ensenmbl.fasta.gz

# Annotation
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz

wget https://github.com/comprna/SUPPA_supplementary_data/raw/master/annotation/Homo_sapiens.GRCh37.75.formatted.gtf.gz
gunzip Homo_sapiens.GRCh37.75.formatted.gtf.gz

wget https://github.com/comprna/SUPPA_supplementary_data/raw/master/annotation/Homo_sapiens.GRCh37.75.formatted.gff.gz
gunzip Homo_sapiens.GRCh37.75.formatted.gff.gz

cd ..

# Samples
mkdir -p Samples
echo "Samples..."
mkdir -p Samples
mv Samples.tar.gz Samples
cd Samples
tar xfz Samples.tar.gz

cd $CF
