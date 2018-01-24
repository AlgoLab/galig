#/bin/bash

ScriptsFold=$(dirname $0)/Scripts

Refs=$1/SimData/Refs
Samples=$1/SimData/Samples/$2
Out=$1/SimEvents/$2/spl_hisat2/Aligns
HisatFold=$1/Tools/hisat2

mkdir -p ${Out}

Index=${Out}/HisatIndex/HI

while read chr
do
    ref=${Refs}/${chr}.fasta
    Sample=${Samples}/${chr}.fasta
    
    out=${Out}/${chr}

    mkdir -p $(dirname ${Index})

    echo "* Indexing $chr - $(date +%r)"
    ${HisatFold}/hisat2-build ${ref} ${Index} &> ${out}.log

    echo "* Aligning $gene ($chr) - $(date +%r)"
    ${HisatFold}/hisat2 -f -x ${Index} -U ${Sample} -S ${out}.sam &>> ${out}.log

    bash ${ScriptsFold}/sam2bam.sh ${out}.sam

    rm -r $(dirname ${Index})
done < ./CHRs
