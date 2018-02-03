#!/bin/bash

WF=$1

for sample in SRR354041 SRR354042 SRR354043
do
    echo "*** Creating reduced annotations (sample ${sample}) - $(date +%r)"
    bash ./RealEvents/cutAnnotations.sh ${WF} ${sample}
    echo ""
    echo "*** Creating truth (sample ${sample}) - $(date +%r)"
    bash ./RealEvents/buildTruth.sh ${WF} ${sample}
    echo ""
    echo "*** Running ASGAL (sample ${sample}) - $(date +%r)"
    bash ./RealEvents/runAsgal.sh ${WF} ${sample}
    echo ""
    echo "*** Analyzing results (sample ${sample}) - $(date +%r)"
    bash ./RealEvents/analyzeAsgal.sh ${WF} ${sample}
    echo ""
done
