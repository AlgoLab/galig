#/bin/bash

ScriptsFold=$(dirname $0)/Scripts

Refs=$1/SimData/Refs
Samples=$1/SimData/Samples/$2
Out=$1/SimEvents/$2/spl_star/Aligns
StarFold=$1/Tools/star/source

mkdir -p ${Out}

Index=${Out}/StarIndex

while read chr
do
    ref=${Refs}/${chr}.fasta
    Sample=${Samples}/${chr}.fasta
    
    out=${Out}/${chr}

    mkdir -p ${Index}

    echo "* Indexing $chr - $(date +%r)"
    ${StarFold}/STAR --runMode genomeGenerate --genomeDir ${Index} --genomeFastaFiles ${ref} &> ${out}.log

    echo "* Aligning $chr - $(date +%r)"
    ${StarFold}/STAR --genomeDir ${Index} --readFilesIn ${Sample} &>> ${out}.log

    sed -i 's/nM:i/NM:i/g' Aligned.out.sam
    mv Aligned.out.sam ${out}.sam
    bash ${ScriptsFold}/sam2bam.sh ${out}.sam

    rm Log.out Log.progress.out Log.final.out SJ.out.tab

    rm -r ${Index}
done < ./CHRs
