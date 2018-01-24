#!/bin/bash

ScriptsFold=$(dirname $0)/Scripts

SimDataFold=$1/SimData

Refs=${SimDataFold}/Refs
Annos=${SimDataFold}/Annotations
OutFold=${SimDataFold}/Refs/Genes/

while read chr
do
    echo -n "${chr} "
    mkdir -p ${OutFold}/${chr}
    python3 ${ScriptsFold}/extractGenes.py ${Refs}/${chr}.fasta ${Annos}/${chr} ${OutFold}/${chr}
done < ./CHRs
echo ""
