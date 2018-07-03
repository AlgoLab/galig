[//]: # (Comment)
# Experiments

This page contains information on how to replicate the experiments
described in _ASGAL_ paper.  To create a fully reproducible data
analysis, we used the
[Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow
manager.

**All the scripts have been tested on ubuntu 18.04.**

## Prerequisites
```bash
sudo apt update
sudo apt-get install build-essential cmake git curl wget unzip samtools python python3 python3-pip \
                     python3-biopython python3-biopython-sql python3-pysam python-scipy python-pysam \
                     python-h5py emboss python-cutadapt libboost1.65-all-dev zlib1g-dev liblzma-dev \
                     libjemalloc1 libjemalloc-dev libghc-bzlib-dev snakemake
pip3 install --user gffutils
```

## Simulated Data

Let us download the data and run the experiments in the following folder:
```bash
SimFold='~/asgal_exp/SimData'
mkdir -p ${SimFold}
```

1. move to the folder containing the snakefile
```bash
cd paper/experiments/SimulatedData
```

2. download the tools:
```bash
bash getTools.sh ${SimFold}
```

3. download the files from [here](https://drive.google.com/open?id=1mbEYLIn9193WdSEBp3rsEpiqC9mS2Vi2) and move them to _SimFold_

4. setup the input files (this script will split the input files and create the reduced annotations):
```bash
bash setupData.sh ${SimFold}
```

5. change the _root_ folder in the ```config.yaml``` file with the desired folder (i.e. ```${SimFold}```)

6. run the experiments using snakemake
```
# check if everything is okay
snakemake -n all
# run the experiments
snakemake all
```

The outputs of the tools are stored in the folder: ``` ${SimFold}/Results ```.

In the same folder, you can find 4 _csv_ which summarize the results:
   * _alignmentsAccuracy.csv_ contains the results on the basewise accuracy of ASGAL aligner and STAR
   * _alignmentsStatistics.csv_ contains the statistics on number of unique primary alignments, number of mismatches and read truncations
   * _full-annot-comparisons.csv_ contains the number of true positives, false positives, and false negatives reported by each tool (ASGAL, SplAdder, rMATS, and SUPPA) when considering the original gene annotations
   * _full-novel-comparisons.csv_ contains the number of true positives, false positives, and false negatives reported by each tool (ASGAL, SplAdder, rMATS, and SUPPA) when considering the reduced annotations


## Real Data

Let us download the data and run the experiments in the following folder:
```bash
RealFold='~/asgal_exp/RealData'
mkdir -p ${RealFold}
```

1. move to the folder containing the snakefile
```bash
cd paper/experiments/RealData
```

2. setup the tools folder (we will create a symbolic link to the _Tools_ folder used before):
```bash
ln -s ${SimFold}/Tools/ ${RealFold}/Tools/
```

3. download the files from [here](https://drive.google.com/open?id=1N5zg3z9XQiOuzpEZUW0HxfKxTGrf5Vtz) and move them to _RealFold_ *(currently, we are uploading the new files)*

4. download the other required data and setup all the data:
```bash
bash setupData.sh ${RealFold}
```

5. change the _root_ folder in the _config.yaml_ file with the desired folder (i.e. ```${RealFold}```)

6. run the experiments using snakemake <!-- ~22779 jobs -->
```
# check if everything is okay
snakemake -n all
# run the experiments
snakemake all
```

The outputs of the tools are stored in the folder: ```${RealFold}/Results ```.