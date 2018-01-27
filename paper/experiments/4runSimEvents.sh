#!/bin/bash

WF=$1

for sample in 5M #10M
do
    echo "*** Creating truth and reduced annotations (sample ${sample}) - $(date +%r)"
    bash ./SimEvents/createTruthAndReducedAnnos.sh ${WF}
    echo ""

    echo "*** Running ASGAL (sample ${sample}) - $(date +%r)"
    bash ./SimEvents/runAsgal.sh ${WF} ${sample}
    echo ""

    echo "*** Running SplAdder with Hisat2 (sample ${sample}) - $(date +%r)"
    echo "** Running Hisat2 (sample ${sample}) - $(date +%r)"
    bash ./SimEvents/runHisat2.sh ${WF} ${sample}
    echo "** Running SplAdder (sample ${sample}) - $(date +%r)"
    bash ./SimEvents/runSpladder.sh ${WF} ${sample} hisat2
    echo ""

    echo "*** Running SplAdder with STAR (sample ${sample}) - $(date +%r)"
    echo "** Running STAR (sample ${sample}) - $(date +%r)"
    bash ./SimEvents/runStar.sh ${WF} ${sample}
    echo "** Running SplAdder (sample ${sample}) - $(date +%r)"
    bash ./SimEvents/runSpladder.sh ${WF} ${sample} star
    echo ""

    echo "*** Analyzing results (sample ${sample}) - $(date +%r)"
    bash ./SimEvents/runAnalyzer.sh ${WF} ${sample}
    echo ""
done
