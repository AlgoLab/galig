#!/bin/bash

ScriptsFold=$(dirname $0)/Scripts

WorkFold=$1/RealData/$2
ResFold=$1/Results/RealEvents/asgal/$2

mkdir -p ${ResFold}

res=${ResFold}/Results.csv
Rres=${ResFold}/ResultsRecall.csv
Pres=${ResFold}/ResultsPrecision.csv
Ires=${ResFold}/IntronsAnalysis.csv
rm -f ${res}
rm -f ${Pres}
rm -f ${Rres}
rm -f ${Ires}

for events in $(ls ${WorkFold}/*.events)
do
    gene=$(basename ${events} .events)
    truth=${WorkFold}/Truth/${gene}.ort.events
    allEvents=${WorkFold}/Truth/${gene}.all.events
    bam=${WorkFold}/${gene}.star.bam
    samtools index ${bam}
    chr=$(head -1 ${WorkFold}/${gene}.gtf | cut -f 1 -d$'\t')
    python3 ${ScriptsFold}/check.py ${events} ${truth} ${allEvents} ${chr} ${gene} >> ${res}
    python3 ${ScriptsFold}/checkRec.py ${events} ${truth} ${chr} ${gene} >> ${Rres}
    python3 ${ScriptsFold}/checkPrec.py ${events} ${truth} ${allEvents} ${chr} ${gene} >> ${Pres}
    python3 ${ScriptsFold}/compareIntrons.py ${events} ${bam} ${chr} ${gene} >> ${Ires}
done
