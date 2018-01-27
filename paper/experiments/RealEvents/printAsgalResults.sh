#!/bin/bash

ResFold=$1

function print {
    sample=$1

    echo "*** ${sample} ***"
    res=${ResFold}/${sample}/Results.csv
    Rres=${ResFold}/${sample}/ResultsRecall.csv
    Ires=${ResFold}/${sample}/IntronsAnalysis.csv

    for ev in ES A3 A5 IR
    do
        totR=$(grep -c $ev $Rres)
        foundR=$(grep $ev $Rres | grep -c " 1$")
        R=$(bc <<< "scale = 3; ($foundR / $totR)")
        echo -e "$ev\t$foundR/$totR=$R"
    done

    echo ""
    FoundRealIntrons=$(cut -f 3 -d',' ${Ires} | paste -sd+ | bc)
    FoundIntrons=$(cut -f 4 -d',' ${Ires} | paste -sd+ | bc)
    RealIntrons=$(cut -f 5 -d',' ${Ires} | paste -sd+ | bc)
    Prec=$(bc <<< "scale = 3; ($FoundRealIntrons / $FoundIntrons)")
    Rec=$(bc <<< "scale = 3; ($FoundRealIntrons / $RealIntrons)")
    
    echo -e "* Introns: $FoundRealIntrons/$FoundIntrons ($Prec)"
    echo ""
}

print SRR354041
print SRR354042
print SRR354043
