#!/bin/bash

if [-s /data/transcripts.txt]
then
  if [-s /data/sample2.fa]
  then
    ./asgal --multi \
        -g /data/genome.fasta \
        -a [annotation.gtf] \
        -s [sample1.fa] \
        -s2 [sample2.fa] \
        -t [transcripts.fasta] \
        -o /data
   else
     ./asgal --multi \
        -g /data/genome.fasta \
        -a [annotation.gtf] \
        -s [sample1.fa] \
        -t [transcripts.fasta] \
        -o /data
   fi
else
  ./asgal -g /data/genome.fasta -a [annotation] -s [sample] -o /data
fi
