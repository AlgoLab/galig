#!/bin/bash

WF=$1

for sample in 5M #10M
do
    echo "*** Creating truth and reduced annotations"
    bash ./SimEvents/createTruthAndReducedAnnos.sh ${WF}
    echo ""

    echo "*** Running ASGAL"
    bash ./SimEvents/runAsgal.sh ${WF} ${sample}
    echo ""

    echo "*** Running SplAdder with Hisat2"
    echo "** Running Hisat2"
    bash ./SimEvents/runHisat2.sh ${WF} ${sample}
    echo "** Running SplAdder"
    bash ./SimEvents/runSpladder.sh ${WF} ${sample} hisat2
    echo ""

    echo "*** Running SplAdder with STAR"
    echo "** Running STAR"
    bash ./SimEvents/runStar.sh ${WF} ${sample}
    echo "** Running SplAdder"
    bash ./SimEvents/runSpladder.sh ${WF} ${sample} star
    echo ""

    echo "*** Analyzing results"
    bash ./SimEvents/runAnalyzer.sh ${WF} ${sample}
    echo ""
done
