#!/bin/bash

ScriptsFold=$(dirname $0)/Scripts

WorkFold=$1/RealData
SpladderFold=$1/Tools/spladder/python

s1=SRR354041
s2=SRR354042
s3=SRR354043

function SingleSpl {
    tool=$1
    log=${WorkFold}/SplAdder/Single/${tool}.log
    rm -f ${log}

    for gtf in $(ls ${WorkFold}/${s1}/*.new.gtf)
    do
        gene=$(basename ${gtf} .new.gtf)
        chr=$(head -1 ${gtf} | cut -f 1 -d$'\t')
        echo "* $gene ($chr) $(date +%r)"

        out=${WorkFold}/SplAdder/Single/${tool}/${gene}
        mkdir -p ${out}
        rm -rf ${out}/spladder

        bam1=${WorkFold}/${s1}/${gene}.${tool}.bam
   	    bam2=${WorkFold}/${s2}/${gene}.${tool}.bam
        bam3=${WorkFold}/${s3}/${gene}.${tool}.bam

	    samtools view -h ${bam1} > sam
	    samtools view ${bam2} >> sam
	    samtools view ${bam3} >> sam
	    samtools view -bS sam > bam
	    samtools sort bam -o bam.sort.bam
	    mv bam.sort.bam ${out}.bam
	    samtools index ${out}.bam
	    rm -f sam bam

	    rm -f ${gtf}.pickle
        python ${SpladderFold}/spladder.py -b ${out}.bam -a ${gtf} -o ${out} &>> ${log}
        cat ${out}/*.txt > ${out}.splout 2> /dev/null
        python3 ${ScriptsFold}/formatSplOut.py ${out}.splout ${gtf} | sort -u > ${out}.events
    done
}

echo "*** Running SplAdder Single + hisat2 $(date +%r)"
SingleSpl hisat2
echo ""
echo "*** Running SplAdder Single + star $(date +%r)"
SingleSpl star
echo ""
