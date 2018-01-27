#!/bin/bash

ScriptsFold=$(dirname $0)/Scripts

ResFold=$1

function printSgal {
    tool=asgal
    echo "*** ${tool} ***"
    echo -e "P\tR\tF"
    python3 ${ScriptsFold}/printPR.py ${ResFold}/${tool}.csv
    
    time=$(python3 ${ScriptsFold}/printTime.py ${ResFold}/${tool}.time)
    ram=$(python3 ${ScriptsFold}/printRAM.py ${ResFold}/${tool}.ram)

    echo -e "\n* Total Time (s):\t\t${time}"
    echo -e "* Max RAM (MB):\t\t\t${ram}"
    echo ""
}

function print {
    tool=$1

    echo "*** ${tool} ***"
    echo -e "P\tR\tF"
    python3 ${ScriptsFold}/printPR.py ${ResFold}/${tool}.csv

    time=$(python3 ${ScriptsFold}/printTime.py ${ResFold}/${tool}.index.time)    
    ram=$(python3 ${ScriptsFold}/printRAM.py ${ResFold}/${tool}.index.ram)
    echo -e "\n* Total Time - Index (s):\t${time}"
    echo -e "* Max RAM - Index (MB):\t\t${ram}"
    
    
    time=$(python3 ${ScriptsFold}/printTime.py ${ResFold}/${tool}.align.time)
    ram=$(python3 ${ScriptsFold}/printRAM.py ${ResFold}/${tool}.align.ram)
    echo -e "\n* Total Time - Align (s):\t${time}"
    echo -e "* Max RAM - Align (MB):\t\t${ram}"
    echo ""
}

printSgal
print hisat2
print star
