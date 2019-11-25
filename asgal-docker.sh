#!/bin/bash

if [ -s /data/transcripts.fa ]
then
    if [ -s /data/sample_2.fa ]
    then
	./asgal --multi \
		-g /data/genome.fa \
		-a /data/annotation.gtf \
		-s /data/sample_1.fa \
		-s2 /data/sample_2.fa \
		-t /data/transcripts.fa \
		-o /data/output
    else
	./asgal --multi \
		-g /data/genome.fasta \
		-a /data/annotation.gtf \
		-s /data/sample_1.fa \
		-t /data/transcripts.fa \
		-o /data/output
    fi
else
    ./asgal -g /data/genome.fasta -a /data/annotation.gtf -s /data/sample_1.fa -o /data/output
fi
