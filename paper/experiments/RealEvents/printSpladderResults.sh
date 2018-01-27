#!/bin/bash

ResFold=$1

function print {
    csv=$1

    echo "*** $(basename ${csv}) ***"

    for ev in ES A3 A5 IR
    do
        totR=$(grep -c $ev $csv)
        foundR=$(grep $ev $csv | grep -c " 1$")
        R=$(bc <<< "scale = 3; ($foundR / $totR)")
        echo -e "$ev\t$foundR/$totR=$R"
    done

    echo ""
}

for csv in $(ls ${ResFold}/*.csv)
do
    print ${csv}
done
