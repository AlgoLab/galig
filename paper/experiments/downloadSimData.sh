#!/bin/bash

CF=$(pwd)

WF=$1
DataFold=${WF}/SimData

mkdir -p ${DataFold}

log="${WF}/getSimData.log"

cd ${DataFold}
# Reference
echo "*** Reference"
echo "** Downloading reference"
#wget http://public.bmi.inf.ethz.ch/projects/2015/spladder/simulation/genome/GRCh37.p13.genome.fa.gz &> ${log}
bash ${CF}/getFromDrive.sh 'https://drive.google.com/open?id=1YbNsL6SXiWYSdVK0oPKuhvO4RoMytY-u' &> ${log}
echo "** Extracting reference"
gunzip -f GRCh37.p13.genome.fa.gz &>> ${log}
echo ""

# Annotation
echo "*** Annotation"
echo "** Downloading annotation"
#wget http://public.bmi.inf.ethz.ch/projects/2015/spladder/annotation/gencode.v19.annotation.mult_iso_subsample_1000_genes.gtf &> ${log}
bash ${CF}/getFromDrive.sh 'https://drive.google.com/open?id=1m33itbejiWmruQ8HUT8j1FOWImwcORDv' &>> ${log}
echo ""

# Samples
echo "*** Samples"
echo "** Downloading sample 5M"
#wget http://public.bmi.inf.ethz.ch/projects/2015/spladder/simulation/reads/5000000_reads.fastq.gz &> ${log}
bash ${CF}/getFromDrive.sh 'https://drive.google.com/open?id=1oef4qtlDJkMTf-K9ueZ3Mw21o8rDpbzS' &>> ${log}
echo "** Extracting sample 5M"
gunzip -f 5000000_reads.fastq.gz &>> ${log}

echo "** Downloading sample 10M"
#wget http://public.bmi.inf.ethz.ch/projects/2015/spladder/simulation/reads/10000000_reads.fastq.gz &> ${log}
bash ${CF}/getFromDrive.sh 'https://drive.google.com/open?id=1xHyKGTeCspZEHqmlGh4OAqWRdW9oUmLu' &>> ${log}
echo "** Extracting sample 10M"
gunzip -f 10000000_reads.fastq.gz &>> ${log}

echo ""

cd ${CF}
