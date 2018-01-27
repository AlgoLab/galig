#!/bin/bash

InFold=$1

ScriptsFold="/pico/scratch/userexternal/ldenti00/SGALScripts"
RefsFold="/gpfs/work/ELIX2_prj16_0/Refs"

#echo -e "\n--- Preparing input ( $(date +%r) ) ---\n"
#bash 1cleanGTFs.sh ${InFold} ${ScriptsFold}

echo -e "\n\n--- Cutting annotations ( $(date +%r) ) ---\n"
bash ${ScriptsFold}/2cutAnnotations.sh ${InFold}

echo -e "\n\n--- Building truth ( $(date +%r) ) ---\n"
bash ${ScriptsFold}/3buildTruth.sh ${InFold} ${ScriptsFold}

echo -e "\n\n--- Running SGAL ( $(date +%r) ) ---\n"
bash ${ScriptsFold}/4runSgal.sh ${InFold} ${RefsFold} ${ScriptsFold}

echo -e "\n\n--- Analyzing ( $(date +%r) ) ---\n"
bash ${ScriptsFold}/5check.sh ${InFold} ${ScriptsFold}
