#!/bin/bash

ScriptsFold=$(dirname $0)/Scripts

Annos=$1/SimData/Annotations
Samples=$1/SimData/Samples/5M
Aligns=$1/AlComp

ResFold=$1/Results/AlComp

mkdir -p ${ResFold}

function analyzeSgal {
    echo "* Analyzing asgal"
    rm -f ${ResFold}/asgal.*

    while read chr
    do
        python3 ${ScriptsFold}/analyzeAligns.py ${Samples}/${chr}.fasta ${Annos}/${chr} ${Aligns}/asgal/${chr} >> ${ResFold}/asgal.csv
    done < ./CHRs

    grep -h "Elapsed (wall clock) time" ${Aligns}/asgal/*/*.time | cut -f 8 -d' ' > ${ResFold}/asgal.time
    grep -h "Maximum resident set size" ${Aligns}/asgal/*/*.time | cut -f 6 -d' ' > ${ResFold}/asgal.ram
}

function analyze {
    tool=$1

    echo "* Analyzing ${tool}"

    rm -f ${ResFold}/${tool}.*
    while read chr
    do
        python3 ${ScriptsFold}/analyzeAligns.py ${Samples}/${chr}.fasta ${Annos}/${chr} ${Aligns}/${tool}/${chr} >> ${ResFold}/${tool}.csv
    done < ./CHRs

    grep -h "Elapsed (wall clock) time" ${Aligns}/${tool}/*/*.index.time | cut -f 8 -d' ' > ${ResFold}/${tool}.index.time
    grep -h "Maximum resident set size" ${Aligns}/${tool}/*/*.index.time | cut -f 6 -d' ' > ${ResFold}/${tool}.index.ram
    grep -h "Elapsed (wall clock) time" ${Aligns}/${tool}/*/*.align.time | cut -f 8 -d' ' > ${ResFold}/${tool}.align.time
    grep -h "Maximum resident set size" ${Aligns}/${tool}/*/*.align.time | cut -f 6 -d' ' > ${ResFold}/${tool}.align.ram
}

analyzeSgal
analyze hisat2
analyze star
