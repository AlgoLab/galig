#!/bin/bash

inFasta=$1
outDir=$2
log=$3

mkdir -p ${outDir}

sample=$(dirname ${inFasta})/$(basename ${inFasta} .fastq)

echo -e "* Converting to FASTA"

cutadapt -o ${sample}.fasta ${sample}.fastq &> ${log}

rm ${sample}.fastq

echo -e "* Trimming sample"
cutadapt -m 40 -a "A{100}" -g "T{100}" ${sample}.fasta > ${sample}.trimmed 2>> ${log}
mv ${sample}.trimmed ${sample}.fasta

echo -e "* Splitting sample"
while read chr
do
    echo -n "${chr} "
    SingleFasta=${outDir}/${chr}.fasta

    # This works since the reads are not splitted on more lines
    grep -A 1 ${chr} ${sample}.fasta > ${SingleFasta}
done < ./CHRs
echo ""
