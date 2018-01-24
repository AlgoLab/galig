#/bin/bash

Refs=$1/SimData/Refs
Annos=$1/SimData/ReducedAnnotations
Samples=$1/SimData/Samples/$2
Out=$1/SimEvents/$2/asgal

while read chr
do
    Ref=${Refs}/${chr}.fasta
    Sample=${Samples}/${chr}.fasta
    ChrAnnos=${Annos}/${chr}
    for f in $(ls -d ${ChrAnnos}/*/)
    do
        gene=$(basename $f)
        for gtf in $(ls $f/*.gtf)
        do
            run=$(basename $gtf .gtf)
            echo "* $gene - $run ($chr) - $(date +%r)"
            out=${Out}/${chr}/${gene}/${run}
            mkdir -p $(dirname ${out})
            ../../bin/SpliceAwareAligner -g ${Ref} -a ${gtf} -s ${Sample} -o ${out}.mem &> ${out}.log
            \time -v -o ${out}.time python3 ../../scripts/detectEvents.py -g ${Ref} -a ${gtf} -m ${out}.mem -o ${out}.events &>> ${out}.log
        done
    done
done < ./CHRs
