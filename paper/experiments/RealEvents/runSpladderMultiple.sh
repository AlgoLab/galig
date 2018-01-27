#!/bin/bash

ScriptsFold=$(dirname $0)/Scripts

WorkFold=$1/RealData
SpladderFold=$1/Tools/spladder/python

s1=SRR354041
s2=SRR354042
s3=SRR354043

function MultSpl {
    tool=$1
    log=${WorkFold}/SplAdder/Multiple/${tool}.log
    rm -f ${log}

    for gtf in $(ls ${WorkFold}/${s1}/*.new.gtf)
    do
        gene=$(basename ${gtf} .new.gtf)
        chr=$(head -1 ${gtf} | cut -f 1 -d$'\t')
        echo "* SplAdder + ${tool} on $gene ($chr) - $(date +%r)"

        out=${WorkFold}/SplAdder/Multiple/${tool}/${gene}
        mkdir -p ${out}
        rm -rf ${out}/spladder

        bam1=${WorkFold}/${s1}/${gene}.${tool}.bam
        bam2=${WorkFold}/${s2}/${gene}.${tool}.bam
        bam3=${WorkFold}/${s3}/${gene}.${tool}.bam

	    rm -f ${gtf}.pickle
        python ${SpladderFold}/spladder.py -b ${bam1},${bam2},${bam3} -a ${gtf} -o ${out} &>> ${log}
        cat ${out}/*.txt > ${out}.splout 2> /dev/null
        python3 ${ScriptsFold}/formatSplOut.py ${out}.splout ${gtf} | sort -u > ${out}.events
    done
}

echo "*** Running SplAdder Multiple + hisat2 $(date +%r)"
MultSpl hisat2
echo ""
echo "*** Running SplAdder Multiple + star $(date +%r)"
MultSpl star
echo ""
