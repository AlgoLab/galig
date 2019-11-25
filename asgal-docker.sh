#!/bin/bash

if [ -s /data/transcripts.fa ]
then
    if [ -s /data/sample_2.fa ]
    then
	/galig/asgal --multi \
		     -g /data/genome.fa \
		     -a /data/annotation.gtf \
		     -s /data/sample_1.fa \
		     -s2 /data/sample_2.fa \
		     -t /data/transcripts.fa \
		     -o /data/output
    else
	/galig/asgal --multi \
		     -g /data/genome.fa \
		     -a /data/annotation.gtf \
		     -s /data/sample_1.fa \
		     -t /data/transcripts.fa \
		     -o /data/output
    fi
else
    /galig/asgal -g /data/genome.fa \
		 -a /data/annotation.gtf \
		 -s /data/sample_1.fa \
		 -o /data/output
fi
