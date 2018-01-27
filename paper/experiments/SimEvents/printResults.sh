#!/bin/bash

ScriptsFold=$(dirname $0)/Scripts

ResFold=$1

for file in $(ls ${ResFold}/*.csv)
do
    fname=$(basename $file .csv)
    echo "*** $fname ***"
    echo -e "Ev\tP\tR\tF"
    for ev in ES A3 A5 IR
    do
        TP=$(grep $ev $file | cut -f 6 -d',' | paste -sd+ | bc)
        FN=$(grep $ev $file | cut -f 7 -d',' | paste -sd+ | bc)
        FP=$(grep $ev $file | cut -f 8 -d',' | paste -sd+ | bc)
        P=$(bc <<< "scale = 3; ($TP / ($TP+$FP))")
        R=$(bc <<< "scale = 3; ($TP / ($TP+$FN))")
        F=$(bc <<< "scale = 3; (2 * ($P * $R) / ($P + $R))")
        echo -e "$ev\t$P\t$R\t$F"
    done
    for ev in "2,-" "3,-" "4,-" "5,-" 
    do
        TP=$(grep $ev $file | cut -f 6 -d',' | paste -sd+ | bc)
        FN=$(grep $ev $file | cut -f 7 -d',' | paste -sd+ | bc)
        FP=$(grep $ev $file | cut -f 8 -d',' | paste -sd+ | bc)
        P=$(bc <<< "scale = 3; ($TP / ($TP+$FP))")
        R=$(bc <<< "scale = 3; ($TP / ($TP+$FN))")
        F=$(bc <<< "scale = 3; (2 * ($P * $R) / ($P + $R))")
        ev=$(echo $ev | cut -f 1 -d',')
        echo -e "$ev\t$P\t$R\t$F"
    done
    echo ""
    time=$(python3 ${ScriptsFold}/printTime.py ${ResFold}/${fname}.time)
    echo -e "* Average Time (s):\t${time}"

    ram=$(python3 ${ScriptsFold}/printRAM.py ${ResFold}/${fname}.ram)
    echo -e "* Average RAM (MB):\t${ram}"
    echo ""
done
