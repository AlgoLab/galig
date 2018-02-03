#!/bin/bash

WF=$1

echo "*** Splitting chromosomes into genes"
bash ./AlComp/extractGenes.sh ${WF}
echo ""

echo "*** Running ASGAL"
bash ./AlComp/runAsgal.sh ${WF}
echo ""

echo "*** Running Hisat2"
bash ./AlComp/runHisat2.sh ${WF}

echo "*** Running STAR"
bash ./AlComp/runStar.sh ${WF}
echo ""

echo "*** Analyzing results"
bash ./AlComp/runAnalyzer.sh ${WF}
echo ""
