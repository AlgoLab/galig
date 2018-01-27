#!/bin/bash

WorkFold=$1/RealData/$2

for gtf in $(ls ${WorkFold}/*.gtf)
do
    gene=$(basename ${gtf} .gtf)
    bed=${WorkFold}/${gene}.bed

    echo -e "* ${gene}"
    trIDs=$(cut -f 4 ${bed} | cut -f 1 -d'-')

    newGTF=${WorkFold}/${gene}.new.gtf
    cat ${gtf} > ${newGTF}
    for id in ${trIDs}
    do
        echo -n "${id} "
        grep -v ${id} ${newGTF} > tmp
        mv tmp ${newGTF}
    done
    if grep -qsP "\ttranscript\t" ${newGTF}
    then
        if grep -qsP "\texon\t" ${newGTF}
        then
            echo ": Ok"
        else
            echo ": 0 exons - The reduced annotation will not be created"
            rm ${newGTF}
        fi
    else
        echo ": 0 transcripts - The reduced annotation will not be created"
        rm ${newGTF}
    fi
    echo ""
done
