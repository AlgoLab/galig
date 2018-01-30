#!/bin/bash

CurrFold=$(pwd)

WF=$1

mkdir -p ${WF}/RealData
cd ${WF}/RealData

log="./RealDataDownload.log"

sample="SRR354041"
echo "* Downloading and extracting ${sample}"
bash ${CurrFold}/getFromDrive.sh 'https://drive.google.com/open?id=1VfAMjnFIcVSdyp-MaXU7FkW4VgLVKHew' &> ${log}
tar xfz ${sample}.tar.gz
rm ${sample}.tar.gz

sample="SRR354042"
echo "* Downloading and extracting ${sample}"
bash ${CurrFold}/getFromDrive.sh 'https://drive.google.com/open?id=1nF-wk9ML1DZvz7PvOrJMueJnl3yycX6c' &>> ${log}
tar xfz ${sample}.tar.gz
rm ${sample}.tar.gz

sample="SRR354043"
echo "* Downloading and extracting ${sample}"
bash ${CurrFold}/getFromDrive.sh 'https://drive.google.com/open?id=1EWKatEx3WaXFsMn3Pr2xTYx8XEnJov4o' &>> ${log}
tar xfz ${sample}.tar.gz
rm ${sample}.tar.gz

cd ${CurrFold}
