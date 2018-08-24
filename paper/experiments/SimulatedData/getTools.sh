#!/bin/bash

CD=$(pwd)

WD=$1

mkdir -p $WD/Tools
cp -r tools/* $WD/Tools

cd $WD/Tools

# ASGAL
git clone --recursive https://github.com/AlgoLab/galig.git
cd galig
make lemon
make sdsl
make
cd ..

# STAR
wget https://github.com/alexdobin/STAR/archive/2.5.4b.tar.gz
tar xfz 2.5.4b.tar.gz
rm 2.5.4b.tar.gz

# SplAdder
git clone https://github.com/ratschlab/spladder.git
cd spladder
git checkout e4ddba022a85cc25d84b37844b1d808a9c1e54b7
cd ..

# rMATS
wget https://sourceforge.net/projects/rnaseq-mats/files/MATS/rMATS.4.0.2.tgz
tar xfz rMATS.4.0.2.tgz
rm rMATS.4.0.2.tgz

# SUPPA2
pip3 install --user SUPPA==2.3

# Salmon
wget https://github.com/COMBINE-lab/salmon/archive/v0.10.2.tar.gz
tar xfz v0.10.2.tar.gz
rm v0.10.2.tar.gz
cd salmon-0.10.2
mkdir build
cd build
cmake -DFETCH_BOOST=TRUE ..
make
make install
cd ../..

# gffread
wget http://ccb.jhu.edu/software/stringtie/dl/gffread-0.9.12.Linux_x86_64.tar.gz
tar xfz gffread-0.9.12.Linux_x86_64.tar.gz
rm gffread-0.9.12.Linux_x86_64.tar.gz

# astalavista
wget http://sammeth.net/data/public/astalavista-4.0.tgz
tar xfz astalavista-4.0.tgz
rm astalavista-4.0.tgz

cd $CD
