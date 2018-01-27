#!/bin/bash

ScriptsFold=$(dirname $0)/Scripts

WorkFold=$1/RealData
ResFold=$1/Results/RealEvents/spladder

mkdir -p ${ResFold}

for fold in $(ls -d ${WorkFold}/SplAdder/*/*/)
do
    run=$(basename ${fold})
    type=$(basename $(dirname ${fold}))
    res=${ResFold}/${run}_${type}.recall.csv
    rm -f ${res}
    for truth in $(ls ${WorkFold}/SRR354041/Truth/*.ort.events)
    do
        gene=$(basename ${truth} .ort.events)
        events=${fold}/${gene}.events
        chr=$(head -1 ${WorkFold}/SRR354041/${gene}.new.gtf | cut -f 1 -d$'\t')
        python3 ${ScriptsFold}/checkRec.py ${events} ${truth} ${chr} ${gene} >> ${res}
    done
done
