#!/bin/bash

WF=$1
mkdir -p ${WF}

# Downloading tools
bash 1downloadTools.sh ${WF}

# Simulated Data
bash 2downloadSimData.sh ${WF}

# Aligners comparison on Simulated Data
bash 3runAlignersComparison.sh ${WF}

# Events detection comparison on Simulated Data
bash 4runSimEvents.sh ${WF}

# Real Data
bash 5downloadRealData.sh ${WF}

# ASGAL on Real Data
bash 6runRealAsgal.sh ${WF}

# SplAdder on Real Data
bash 7runRealSplAdder.sh ${WF}
