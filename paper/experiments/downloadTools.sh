#!/bin/bash

RED='\033[0;31m'
NC='\033[0m'

RF=$(pwd)
ToolsFold=$1/Tools
log="./prereqs.log"

function installSTAR {
    echo -e "\n* STAR"
    echo -e "********"

    rm -rf ./Star

    echo -e "* Cloning..."
    git clone https://github.com/alexdobin/STAR.git &>> ${log}
    mv STAR star

    echo -e "* Building..."
    cd ./star/source &>> ${log}
    if make STAR &>> ${log} ; then
        echo -e "* STAR builds correctly"
    else
        echo -e "${RED}* Error in building STAR${NC}"
        exit 1
    fi

    echo -e "\n\n" >> ${log}

    cd ${ToolsFold}
}

function installHISAT2 {
    echo -e "\n* Hisat2"
    echo -e "**********"
    hisatF="hisat2-2.1.0"

    rm -rf ./Hisat2
    echo -e "* Downloading..."
    wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/${hisatF}-source.zip &>> ${log}

    echo -e "* Unzipping..."
    unzip ${hisatF}-source.zip &>> ${log}
    rm ${hisatF}-source.zip &>> ${log}
    mv ${hisatF} ./hisat2 &>> ${log}
    cd ./hisat2 &>> ${log}

    echo -e "* Building..."
    if make &>> ${log} ; then
        echo -e "* Hisat2 builds correctly"
    else
        echo -e "${RED}* Error in building Hisat2${NC}"
        exit 1
    fi

    echo -e "\n\n" >> ${log}

    cd ${ToolsFold}
}

function installSPL {
    echo -e "\n* SplAdder"
    echo -e "************"

    rm -rf ./Spladder

    echo -e "* Cloning..."
    git clone https://github.com/ratschlab/spladder.git &>> ${log}

    chmod +x spladder/python/spladder.py

    echo -e "\n\n" >> ${log}

    cd ${ToolsFold}
}

function installASTA {
    echo -e "\n* Astalavista"
    echo -e "***************"

    rm -rf ./Astalavista

    echo -e "* Downloading..."
    wget http://sammeth.net/data/public/astalavista-4.0.tgz &>> ${log}
    tar xvfz astalavista-4.0.tgz &>> ${log}
    rm astalavista-4.0.tgz
    mv astalavista-4.0 astalavista &>> ${log}

    echo -e "\n\n" >> ${log}

    cd ${ToolsFold}
}

echo -e "\n\t*******************"
echo -e "\t* Tools Installer *"
echo -e "\t*******************"

mkdir -p ${ToolsFold}

cd ${ToolsFold}

installSTAR
installHISAT2
installSPL
installASTA

cd ${RF}
