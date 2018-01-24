#!/bin/bash

sam=$1

F=$(dirname ${sam})
n=$(basename ${sam} .sam)

samtools view -bS -F 4 ${F}/${n}.sam > ${F}/${n}.bam
samtools sort ${F}/${n}.bam -o ${F}/${n}.sorted.bam
mv ${F}/${n}.sorted.bam ${F}/${n}.bam
samtools index ${F}/${n}.bam
