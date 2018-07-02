#!/bin/bash

WF=$1
log="${WF}/setupSimData.log"
rm -f ${log}

# Reference
echo "*** Splitting reference"
mkdir -p ${WF}/Refs
gunzip ${WF}/GRCh37.p13.genome.fa.gz
python3 ./scripts/splitReference.py ${WF}/GRCh37.p13.genome.fa ${WF}/Refs 2>> ${log}
rm ${WF}/GRCh37.p13.genome.fa
echo ""

# Annotation
echo "*** Splitting annotation"
mkdir -p ${WF}/RealAnnotations
python3 ./scripts/splitAnnotation.py ${WF}/gencode.v19.annotation.mult_iso_subsample_1000_genes.gtf ${WF}/RealAnnotations/ 2>> ${log}
rm -f ${WF}/gencode.v19.annotation.mult_iso_subsample_1000_genes.gtf*
echo ""

# Samples
echo "*** Cleaning sample 5M"
gunzip ${WF}/5000000_reads.fastq
bash ./scripts/cleanAndSplitSample.sh ${WF}/5000000_reads.fastq ${WF}/Samples/5M ${log}
rm -f ${WF}/5000000_reads.fasta

echo "*** Cleaning sample 10M"
gunzip ${WF}/10000000_reads.fastq.gz
bash ./scripts/cleanAndSplitSample.sh ${WF}/10000000_reads.fastq ${WF}/Samples/10M ${log}
rm -f ${WF}/10000000_reads.fasta
echo ""

echo "*** Creating truth and reduced annotations"
Annos=${WF}/RealAnnotations
NewAnnos=${WF}/ReducedAnnotations
ToolsFold=${WF}/Tools/astalavista-4.0/bin

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
do
    for gtf in $(ls ${Annos}/${chr}/*.gtf)
    do
        gene=$(basename ${gtf} .gtf)
        echo -e "\n** ${gene} (${chr}) - $(date +%r)"

        #Creating the reduced annotations...
        echo "* Running ASTALAVISTA"
        out=${NewAnnos}/${chr}/${gene}
        mkdir -p ${out}

        ${ToolsFold}/astalavista -t asta -i ${gtf} -o ${out}.asta.gz &>> ${out}.log
        gunzip -f ${out}.asta.gz

        echo "* Formatting events"
        python3 ./scripts/formatAsta.py ${out}.asta > ${out}.events

        sort -u ${out}.events -o ${out}.events
        python3 ./scripts/addUniqueID.py ${out}.events > tmp
        mv tmp ${out}.events

        echo "* Creating reduced annotations"
        python3 ./scripts/transcriptsRemover.py ${gtf} ${out}.events ${NewAnnos}

        #Combining reduced annotations
        echo "* Combining reduced annotations"
        if [ -n "$(ls -A ${out})" ]
        then
            n=$(ls ${out} | wc -l)
            equalsInfo=${out}/equals.info
            rm -f ${equalsInfo}

            if [ $n = 1 ]
            then
                echo -n "Only one annotation: "
                echo -n "writing truth, "
                ev=$(ls ${out})
                echo E1 $(basename ${ev} .gtf) > ${equalsInfo}
                echo "moving"
                mv ${out}/${ev} ${out}/E1.gtf
            fi
            if (( $n > 1 ));
            then
                echo -n "More than one annotation: "
                echo -n "merging, "
                for f1 in $(ls ${out}/*.gtf)
                do
                    fname1=$(basename $f1 .gtf)
                    for f2 in $(ls ${out}/*.gtf)
                    do
                        fname2=$(basename $f2 .gtf)
                        min=$(echo -e "$fname1\n$fname2" | sort | head -1)
                        max=$(echo -e "$fname1\n$fname2" | sort | tail -1)
                        if ! [[ $min = $max ]]
                        then
                            d=$(diff $f1 $f2 | wc -l)
                            echo $min $max $d >> ${equalsInfo}
                        fi
                    done
                done
                sort -u ${equalsInfo} -o ${equalsInfo}
                python3 ./scripts/combineSame.py ${equalsInfo} > tmp
                mv tmp ${equalsInfo}
                echo "moving."
                python3 ./scripts/moveAnnos.py ${equalsInfo} ${out}
            fi
        else
            echo "Zero events: deleting gene."
            rm -r ${out}
            rm -r ${out}.events
        fi
    done
done
