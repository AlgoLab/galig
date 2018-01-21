[Main Page](index)

# Experiments

--- Work in Progress ---

Before starting, we have to set three variables: the folder where data for simulated experiments will be downloaded, the folder for the output of simulated experiments and a folder for real experiments.
```
ToolsFold=/home/user/tmp/SGALExps/Tools
SimDataFold=/home/user/tmp/SGALExps/SimData
SimOutFold=/home/user/tmp/SGALExps/Out
RealWorkFold=/home/user/tmp/SGALExps/RealData
```

### Prerequisites
1. Update and install prerequisites (apt)
    ```
    sudo apt-get update
    sudo apt-get install build-essential g++ python python3 python3-pip python3-biopython python3-biopython-sql wget git cmake unzip samtools python-pysam python3-pysam zlib1g-dev emboss cutadapt python-h5py python-scipy
    ```
2. Update and install prerequisites (pip3)
    ```
    pip3 install --user gffutils
    ```
3. Install tools used throughout the experiments (ASGAl, SplAdder, Hisat2, STAR, AStalavista) 
    ```
    bash downloadTools.sh ${ToolsFold}
    ```
4. Add tools to path variable
    ```
    PATH=${PATH}:${ToolsFold}/asgal/bin:${ToolsFold}/asgal/scripts
    PATH=${PATH}:${ToolsFold}/hisat2
    PATH=${PATH}:${ToolsFold}/star/source
    PATH=${PATH}:${ToolsFold}/spladder/python
    PATH=${PATH}:${ToolsFold}/astalavista
    ```

### Simulated data
In this analysis...
##### Download data
We will now download...
1. download and unzip the reference, and create a reference for each chromosome
    ```
    wget http://public.bmi.inf.ethz.ch/projects/2015/spladder/simulation/genome/GRCh37.p13.genome.fa.gz
    gunzip -f GRCh37.p13.genome.fa.gz
    python3 DownloadData/splitReference.py GRCh37.p13.genome.fa ${SimDataFold}/Refs
    rm GRCh37.p13.genome.fa
    ```

2. download gene annotations and create one annotation for each gene
    ```
    wget http://public.bmi.inf.ethz.ch/projects/2015/spladder/annotation/gencode.v19.annotation.mult_iso_subsample_1000_genes.gtf
    python3 DownloadData/splitAnnotation.py gencode.v19.annotation.mult_iso_subsample_1000_genes.gtf ${SimDataFold}/Annotations/
    rm gencode.v19.annotation.mult_iso_subsample_1000_genes.gtf*
    ```
3. download RNA-Seq samples, removes poly(A) tails, and divide them into 24 smaller samples, one for each chromosome
    ```
    wget http://public.bmi.inf.ethz.ch/projects/2015/spladder/simulation/reads/5000000_reads.fastq.gz
    gunzip -f 5000000_reads.fastq.gz
    bash DownloadData/cleanAndSplitSample.sh 5000000_reads.fastq ${SimDataFold}/Samples/5M
    rm 5000000_reads.fasta
    
    wget http://public.bmi.inf.ethz.ch/projects/2015/spladder/simulation/reads/10000000_reads.fastq.gz
    gunzip -f 10000000_reads.fastq.gz
    bash DownloadData/cleanAndSplitSample.sh 10000000_reads.fastq ${SimDataFold}/Samples/10M
    rm 10000000_reads.fasta
    ```

##### Aligners comparison
In this analysis...
1. create a reference for each gene
    ```
    bash AlComp/extractGenes.sh ${SimDataFold}/Refs ${SimDataFold}/Annotations/ ${SimDataFold}/Refs/Genes/
    ```
2. run SGAL
    ```
    bash AlComp/runSgal.sh ${SimDataFold}/Refs ${SimDataFold}/Annotations ${SimDataFold}/Samples/5M ${SimOutFold}/AlComp
    ```
3. run STAR and HISAT2
    ```
    bash AlComp/runHisat2.sh ${SimDataFold}/Refs/Genes ${SimDataFold}/Samples/5M ${SimOutFold}/AlComp
    bash AlComp/runStar.sh ${SimDataFold}/Refs/Genes ${SimDataFold}/Samples/5M ${SimOutFold}/AlComp
    ```
4. analyze and print the results (RESULTS FORMATTING)
    ```
    bash AlComp/runAnalyzer.sh ${SimDataFold}/Annotations ${SimDataFold}/Samples/5M ${SimOutFold} Results1
    bash AlComp/printResults.sh Results1
    ```

##### Events detection
In this analysis...
```
SampleName=5M
```
1. create the ground truth using AStalavista and the reduced annotations
    ```
    bash SimEvents/createTruthAndReducedAnnos.sh ${SimDataFold}/Annotations ${SimDataFold}/ReducedAnnotations
    ```
2. run ASGAL
    ```
    bash SimEvents/runSgal.sh ${SimDataFold}/Refs ${SimDataFold}/ReducedAnnotations ${SimDataFold}/Samples/${SampleName} ${SimOutFold}/SimEvents/${SampleName}
    ```
3. run STAR and HISAT2
    ```
    bash SimEvents/runHisat2.sh ${SimDataFold}/Refs ${SimDataFold}/Samples/${SampleName} ${SimOutFold}/SimEvents/${SampleName}
    bash SimEvents/runStar.sh ${SimDataFold}/Refs ${SimDataFold}/Samples/${SampleName} ${SimOutFold}/SimEvents/${SampleName}
    ```
4. run SplAdder
    ```
    bash SimEvents/runSpladder.sh ${SimDataFold}/ReducedAnnotations ${SimOutFold}/SimEvents/${SampleName}
    ```
5. analyze and print results (RESULTS FORMATTING)
    ```
    bash SimEvents/runAnalyzer.sh ${SimDataFold}/ReducedAnnotations ${SimOutFold}/SimEvents/${SampleName} Results2
    bash SimEvents/printResults.sh Results2
    ```

### Real Data
In this analysis...
##### Download data
We will now download...
```
DS=SimpleRealDataset
```
1. download data
    ```
    wget ...
    tar xfz ...
    mv ... ${RealWorkFold}
    rm ...
    ```
##### Events detection
In this analysis...
1. Create the reduced annotations using the differential expressed transcripts and build the truth using AStalavista
    ```
    bash RealEvents/cutAnnotations.sh ${RealWorkFold}/${DS}/
    bash RealEvents/buildTruth.sh ${RealWorkFold}/${DS}/
    ```
2. Run SGAL
    ```
    bash RealEvents/runSgal.sh ${RealWorkFold}/${DS} ${SimDataFold}/Refs
    ```
3. analyze and print results
    ```
    bash RealEvents/runAnalyzer.sh ${RealWorkFold}/${DS} Results3_${DS}
    bash RealEvents/printResults.sh Results3_${DS}
    ```
### Events detection (SplAdder)
In this analysis, ... we need all the sample.
1. run SplAdder
    ```
    bash RealEvents/runSpladder.sh ${RealWorkFold} 3
    bash RealEvents/runSpladder.sh ${RealWorkFold} 0
    bash RealEvents/runSpladderMultiple.sh ${RealWorkFold}
    bash RealEvents/runSpladderCombined.sh ${RealWorkFold}
    ```
2. analyze and print SplAdder results
    ```
    bash RealEvents/analyzeSpladder.sh ${RealWorkFold} RealResultsSplAdder
    bash RealEvents/printSplResults.sh RealResultsSplAdder
    ```
