#!/bin/bash

ScriptsFold=$(dirname $0)/Scripts

sample=$2

NewAnnos=$1/SimData/ReducedAnnotations
Out=$1/SimEvents/${sample}
ResFold=$1/Results/SimEvents

mkdir -p ${ResFold}

# Checking ASGAL
Res=${ResFold}/asgal.${sample}
rm -f ${Res}.csv

while read chr
do
    echo "* Checking $chr (asgal)"
    Truths=${NewAnnos}/${chr}
    Events=${Out}/asgal/${chr}

    for fold in $(ls -d ${Events}/*)
    do
        gene=$(basename ${fold})
        for ev in $(ls ${Events}/${gene}/*.events)
        do
    	    run=$(basename ${ev} .events)
            gtf=${NewAnnos}/${chr}/${gene}/${run}.gtf
            python3 ${ScriptsFold}/check.py ${ev} ${Truths}/${gene}/equals.info ${Truths}/${gene}.events ${gtf} >> ${Res}.csv
        done
    done

    grep -h "Elapsed (wall clock) time" ${Out}/asgal/*/*/*.time | cut -f 8 -d' ' > ${Res}.time
    grep -h "Maximum resident set size" ${Out}/asgal/*/*/*.time | cut -f 6 -d' ' > ${Res}.ram
done < ./CHRs

#Checking SplAdder
function spladder {
    tool=$1
    
    Res=${ResFold}/spladder.${tool}.${sample}
    rm -f ${Res}.csv
    while read chr
    do
        echo "* Checking $chr (spladder + ${tool})"

        Truths=${NewAnnos}/${chr}
        Events=${Out}/spl_${tool}/Events/${chr}

        for fold in $(ls -d ${Events}/*)
        do
            gene=$(basename ${fold})
            for ev in $(ls ${Events}/${gene}/*.events)
            do
                run=$(basename ${ev} .events)
                gtf=${NewAnnos}/${chr}/${gene}/${run}.gtf
                python3 ${ScriptsFold}/check.py ${ev} ${Truths}/${gene}/equals.info ${Truths}/${gene}.events ${gtf} 1 >> ${Res}.csv
            done
        done

        grep -h "Elapsed (wall clock) time" ${Out}/spl_${tool}/Events/*/*/*.time | cut -f 8 -d' ' > ${Res}.time
        grep -h "Maximum resident set size" ${Out}/spl_${tool}/Events/*/*/*.time | cut -f 6 -d' ' > ${Res}.ram
    done < ./CHRs
}

spladder hisat2
spladder star
