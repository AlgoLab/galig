#!/bin/bash

ScriptsFold=$(dirname $0)/Scripts

WorkFold=$1/RealData
SpladderFold=$1/Tools/spladder/python
conf=$2

function runSpladder {
    ds=$1
    tool=$2

    log=${WorkFold}/SplAdder/Conf${conf}/${ds}_${tool}.log
    rm -f ${log}

    for gtf in $(ls ${WorkFold}/${ds}/*.new.gtf)
    do
        gene=$(basename ${gtf} .new.gtf)
        chr=$(head -1 ${gtf} | cut -f 1 -d$'\t')
        echo "* SplAdder${conf} + ${tool} on $gene ($chr) - $(date +%r)"
        bam=${WorkFold}/${ds}/${gene}.${tool}.bam

        out=${WorkFold}/SplAdder/Conf${conf}/${ds}_${tool}/${gene}
        mkdir -p ${out}

        cp ${bam} ${bam}.old
        samtools view -h ${bam} > sam
        sed -i 's/nM:i:/NM:i:/g' sam
        samtools view -bS sam > bam
        samtools sort bam -o bam.sort.bam
        mv bam.sort.bam ${bam}
        samtools index ${bam}
        rm -f sam bam

        rm -f ${gtf}.pickle
        python ${SpladderFold}/spladder.py -b ${bam} -a ${gtf} -o ${out} -c ${conf} &>> ${log}
        cat ${out}/*.txt > ${out}.splout 2> /dev/null
        python3 ${ScriptsFold}/formatSplOut.py ${out}.splout ${gtf} | sort -u > ${out}.events
    done
}

for sample in SRR354041 SRR354042 SRR354043
do
    echo "*** Running SplAdder${conf} + hisat2 - $sample $(date +%r)"
    runSpladder $sample hisat2
    echo ""
    echo "*** Running SplAdder${conf} + star - $sample $(date +%r)"
    runSpladder $sample star
    echo ""
done
