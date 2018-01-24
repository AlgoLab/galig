#/bin/bash

ScriptsFold=$(dirname $0)/Scripts

tool=$3

NewAnnos=$1/SimData/ReducedAnnotations
Out=$1/SimEvents/$2/spl_${tool}/
SpladderFold=$1/Tools/spladder/python

while read chr
do
    Aligns=${Out}/Aligns/${chr}.bam

    for f in $(ls -d ${NewAnnos}/${chr}/*/)
    do
        gene=$(basename $f)
        for gtf in $(ls $f/*.gtf)
        do
            run=$(basename $gtf .gtf)
            echo "* $gene - $run ($chr) - $(date +%r)"
            rm -f ${gtf}.pickle
            out=${Out}/Events/${chr}/${gene}/${run}
            mkdir -p $(dirname ${out})
            rm -fr ${out}/spladder
            \time -v -o ${out}.time python ${SpladderFold}/spladder.py -b ${Aligns} -o ${out} -a ${gtf} &> ${out}.log
            cat ${out}/*.txt > ${out}.splout
            python3 ${ScriptsFold}/formatSplOut.py ${out}.splout ${gtf} | sort -u > ${out}.events
        done
    done
done < ./CHRs
