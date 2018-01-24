#!/bin/bash

ScriptsFold=$(dirname $0)/Scripts

ResFold=$1

function printSgal {
    tool=asgal

    python3 ${ScriptsFold}/printPR.py ${ResFold}/${tool}.csv

    python3 ${ScriptsFold}/printTime.py ${ResFold}/${tool}.time

    python3 ${ScriptsFold}/printRAM.py ${ResFold}/${tool}.ram
}

function print {
    tool=$1

    python3 ${ScriptsFold}/printPR.py ${ResFold}/${tool}.csv

    python3 ${ScriptsFold}/printTime.py ${ResFold}/${tool}.index.time
    python3 ${ScriptsFold}/printTime.py ${ResFold}/${tool}.align.time

    python3 ${ScriptsFold}/printRAM.py ${ResFold}/${tool}.index.ram
    python3 ${ScriptsFold}/printRAM.py ${ResFold}/${tool}.align.ram
}

printSgal
print hisat2
print star
