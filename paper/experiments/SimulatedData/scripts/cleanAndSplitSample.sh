#!/bin/bash

inFasta=$1
outDir=$2
log=$3

mkdir -p ${outDir}

sample=$(dirname ${inFasta})/$(basename ${inFasta} .fastq)

echo -e "* Converting to FASTA"
fastq_to_fasta -i ${sample}.fastq -o ${sample}.fasta

rm ${sample}.fastq

echo -e "* Trimming sample"
cutadapt -m 40 -a "A{100}" -g "T{100}" ${sample}.fasta > ${sample}.trimmed 2>> ${log}
mv ${sample}.trimmed ${sample}.fasta

echo -e "* Splitting sample"
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
    echo -n "chr${chr} "
    SingleFasta=${outDir}/chr${chr}.fasta

    # This works since the reads are not splitted on more lines
    grep -A 1 "^>chr${chr}:" ${sample}.fasta > ${SingleFasta}
done
echo ""
