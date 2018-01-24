#!/bin/bash

ScriptsFold=$(dirname $0)/Scripts

Annos=$1/SimData/Annotations
NewAnnos=$1/SimData/ReducedAnnotations
ToolsFold=$1/Tools/astalavista/bin

while read chr
do
    for gtf in $(ls ${Annos}/${chr}/*.gtf)
    do
        gene=$(basename ${gtf} .gtf)
        echo -e "\n** ${gene} (${chr}) - $(date +%r)"

        #Creating the fake annotations...
        echo "* Running ASTALAVISTA"
        out=${NewAnnos}/${chr}/${gene}
        mkdir -p ${out}

        ${ToolsFold}/astalavista -t asta -i ${gtf} -o ${out}.asta.gz &>> ${out}.log
        gunzip -f ${out}.asta.gz

        echo "* Formatting events"
        python3 ${ScriptsFold}/formatAsta.py ${out}.asta > ${out}.events

        sort -u ${out}.events -o ${out}.events
        python3 ${ScriptsFold}/addUniqueID.py ${out}.events > tmp
        mv tmp ${out}.events

        echo "* Creating fake annotations"
        python3 ${ScriptsFold}/transcriptsRemover.py ${gtf} ${out}.events ${NewAnnos}

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
                python3 ${ScriptsFold}/combineSame.py ${equalsInfo} > tmp
                mv tmp ${equalsInfo}
                echo "moving."
                python3 ${ScriptsFold}/moveAnnos.py ${equalsInfo} ${out}
            fi
        else
            echo "Zero events: deleting gene."
            rm -r ${out}
            rm -r ${out}.events
        fi
    done
done < ./CHRs
