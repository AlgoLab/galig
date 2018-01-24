#/bin/bash

SimDataFold=$1/SimData
OutFold=$1/AlComp/star
ToolsFold=$1/Tools/star/source

Refs=${SimDataFold}/Refs/Genes
Annos=${SimDataFold}/Annotations
Samples=${SimDataFold}/Samples/5M/

mkdir -p ${OutFold}

Index=${OutFold}/index

while read chr
do
    for ref in $(ls ${Refs}/${chr}/*.fasta)
    do
        gene=$(basename ${ref} .fasta)
        chr=$(basename $(dirname ${ref}))

        Sample=${Samples}/${chr}.fasta
    
        out=${OutFold}/${chr}/${gene}
        mkdir -p $(dirname ${out})

        mkdir -p ${Index}
        echo "* Indexing $gene ($chr) - $(date +%r)"
        L=`infoseq -auto -nocolumns -only -noheading -name -length ${ref} | cut -d '|' -f 2`
        NB=`echo ${L} | awk '{printf "%11.9f\n", (log($1)/log(2))/2 -1}' | cut -d '.' -f 1`
        \time -v -o ${out}.index.time ${ToolsFold}/STAR --runMode genomeGenerate --genomeDir ${Index} --genomeFastaFiles ${ref} --genomeSAindexNbases ${NB} &> ${out}.log

        echo "* Aligning $gene ($chr) - $(date +%r)"
        \time -v -o ${out}.align.time ${ToolsFold}/STAR --genomeDir ${Index} --readFilesIn ${Sample} &>> ${out}.log

        sed -i 's/nM:i/NM:i/g' Aligned.out.sam
        mv Aligned.out.sam ${out}.sam

        rm Log.out Log.progress.out Log.final.out SJ.out.tab

        rm -r ${Index}
        echo ""
    done
done < ./CHRs
