#/bin/bash

SimDataFold=$1/SimData
OutFold=$1/AlComp/hisat2
ToolsFold=$1/Tools/hisat2

Refs=${SimDataFold}/Refs/Genes
Annos=${SimDataFold}/Annotations
Samples=${SimDataFold}/Samples/5M/

mkdir -p ${OutFold}

Index=${OutFold}/index/HI

while read chr
do
    for ref in $(ls ${Refs}/${chr}/*.fasta)
    do
        gene=$(basename ${ref} .fasta)
        chr=$(basename $(dirname ${ref}))

        Sample=${Samples}/${chr}.fasta
    
        out=${OutFold}/${chr}/${gene}
        mkdir -p $(dirname ${out})

        mkdir -p $(dirname ${Index})
        echo "* Indexing $gene ($chr) - $(date +%r)"
        \time -v -o ${out}.index.time ${ToolsFold}/hisat2-build ${ref} ${Index} &> ${out}.log

        echo "* Aligning $gene ($chr) - $(date +%r)"
        \time -v -o ${out}.align.time ${ToolsFold}/hisat2 -f -x ${Index} -U ${Sample} -S ${out}.sam &>> ${out}.log
        samtools view -S -F 4 ${out}.sam > ${out}.sam.tmp
        mv ${out}.sam.tmp ${out}.sam
        rm -r $(dirname ${Index})
        echo ""
    done
done < ./CHRs
