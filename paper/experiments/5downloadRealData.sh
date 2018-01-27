#!/bin/bash

CurrFold=$(pwd)

WF=$1

cd ${WF}
mkdir -p RealData
cd RealData

for sample in SRR354041 SRR354042 SRR354043
do
    echo "* Downloading and extracting ${sample}"
    #scp lxc_gralign:/data/SgalRealExp/${sample}.tar.xz .
    #tar xJf ${sample}.tar.xz
    #rm ${sample}.tar.xz
done

cd ${CurrFold}
