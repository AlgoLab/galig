#!/bin/bash

ScriptsFold=$(dirname $0)/Scripts

WorkFold=$1/RealData/$2

AstaFold=$1/Tools/astalavista/bin

#Create ASTA from all gtfs
rm -f ${WorkFold}/AllGenes.gtf
for gtf in $(ls ${WorkFold}/*.gtf | grep -v "new")
do
    cat ${gtf} >> ${WorkFold}/AllGenes.gtf
done

${AstaFold}/astalavista -t asta -i ${WorkFold}/AllGenes.gtf -o ${WorkFold}/AllGenes.asta.gz
gunzip -f ${WorkFold}/AllGenes.asta.gz

fullAsta=${WorkFold}/AllGenes.asta

mkdir -p ${WorkFold}/Truth/

for bed in $(ls ${WorkFold}/*.bed)
do
    gene=$(basename ${bed} .bed)
    out=${WorkFold}/Truth/${gene}

    # All events in the gene
    trIDs=$(grep -P "\ttranscript\t" ${WorkFold}/${gene}.gtf | \
                   egrep -o "transcript_id \"ENST[0-9]*\.[0-9]*" | \
                   cut -f 2 -d"\"")
    rm -f ${out}.all.asta
    for id in ${trIDs}
    do
        grep ${id} ${fullAsta} >> ${out}.all.asta
    done
    python3 ${ScriptsFold}/formatAsta.py ${out}.all.asta > ${out}.all.events.full
    python3 ${ScriptsFold}/filterAllEvents.py ${out}.all.events.full ${WorkFold}/${gene}.gtf | sort -u > ${out}.all.events

    # Events involving only removed transcript
    trIDs=$(cut -f 4 ${bed} | cut -f 1 -d'-')
    rm -f ${out}.ort.asta
    for id in ${trIDs}
    do
        grep ${id} ${fullAsta} >> ${out}.ort.asta
    done
    python3 ${ScriptsFold}/formatAsta.py ${out}.ort.asta > ${out}.ort.events.full
    python3 ${ScriptsFold}/filterEvents.py ${out}.ort.events.full ${WorkFold}/${gene}.gtf ${bed} | sort -u > ${out}.ort.events
done
