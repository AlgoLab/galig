#/bin/bash

SimDataFold=$1/SimData
OutFold=$1/AlComp/asgal

Refs=${SimDataFold}/Refs
Annos=${SimDataFold}/Annotations
Samples=${SimDataFold}/Samples/5M/

mkdir -p ${OutFold}

while read chr
do
    Ref=${Refs}/${chr}.fasta
    Sample=${Samples}/${chr}.fasta
    ChrAnnos=${Annos}/${chr}
    ChrOut=${OutFold}/${chr}
    mkdir -p ${ChrOut}

    for gtf in $(ls ${ChrAnnos}/*.gtf)
    do
        gene=$(basename $gtf .gtf)
        echo "* Aligning $gene ($chr) - $(date +%r)"
        out=${ChrOut}/${gene}
        
        \time -v -o ${out}.time sh -c "../../bin/SpliceAwareAligner -g ${Ref} -a ${gtf} -s ${Sample} -o ${out}.mem ; 
        python3 ../../scripts/formatSAM.py -m ${out}.mem -a ${gtf} -g ${Ref} -o ${out}.sam" &>> ${out}.log
    done
done < ./CHRs
