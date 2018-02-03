#!/bin/bash

WF=$1
mkdir -p ${WF}

# Downloading tools
bash downloadTools.sh ${WF}

# Simulated Data
bash downloadSimData.sh ${WF}
bash setupSimData.sh ${WF}

# Aligners comparison on Simulated Data
bash runAlignersComparison.sh ${WF}

# Events detection comparison on Simulated Data
bash runSimEvents.sh ${WF}

# Real Data
bash downloadRealData.sh ${WF}

# ASGAL on Real Data
bash runRealAsgal.sh ${WF}

# SplAdder on Real Data
bash runRealSplAdder.sh ${WF}
