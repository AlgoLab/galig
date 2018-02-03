#!/bin/bash

WF=$1

# Running SplAdder with hisat2 and STAR with default parameters
bash ./RealEvents/runSpladder.sh ${WF} 3

# Running SplAdder with hisat2 and STAR with 0 confidence
bash ./RealEvents/runSpladder.sh ${WF} 0

# Running SplAdder with hisat2 and STAR, multiple sample
bash ./RealEvents/runSpladderMultiple.sh ${WF}

# Running SplAdder with hisat2 and STAR, single big sample
bash ./RealEvents/runSpladderSingle.sh ${WF}

# Analyzing SplAdder
bash ./RealEvents/analyzeSplAdder.sh ${WF}
