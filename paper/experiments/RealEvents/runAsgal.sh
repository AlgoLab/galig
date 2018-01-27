#!/bin/bash

ScriptsFold=$(dirname $0)/Scripts

WorkFold=$1/RealData/$2
RefsFold=$1/SimData/Refs

for gtf in $(ls ${WorkFold}/*.new.gtf)
do
    gene=$(basename ${gtf} .new.gtf)
    chr=$(head -1 ${gtf} | cut -f 1 -d$'\t')
    echo "* $gene ($chr) $(date +%r)"

    Ref=${RefsFold}/${chr}.fasta
    Sample=${WorkFold}/${gene}.fa

    OutFile=${WorkFold}/${gene}

    ../../bin/SpliceAwareAligner -g ${Ref} -a ${gtf} -s ${Sample} -o ${OutFile}.mem &> ${OutFile}.log
    python3 ../../scripts/detectEvents.py --allevents -g ${Ref} -a ${gtf} -m ${OutFile}.mem -o ${OutFile}.events &>> ${OutFile}.log
done
