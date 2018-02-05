[Main Page](index)

# Experiments

This page contains information on how to replicate the experiments described in _ASGAL_ paper. **All the scripts have been tested only on ubuntu 17.10.**

##### Data availability
All the data used in the experiments are available at the following locations:
* [simulated data](https://drive.google.com/open?id=1mbEYLIn9193WdSEBp3rsEpiqC9mS2Vi2), this is a copy of the data available [here](https://public.bmi.inf.ethz.ch/projects/2015/spladder/)
* [real data](https://drive.google.com/open?id=1N5zg3z9XQiOuzpEZUW0HxfKxTGrf5Vtz)

Notice that our scripts will download and set up the data
automatically. The data are hosted on _google drive_ and to
download them, we use the script available
[here](https://gist.github.com/darencard/079246e43e3c4b97e373873c6c9a3798).

However, it is not guaranteed that this script will work indefinitely: if it fails, please let
me known. A possible workaround is the following one:

1. download manually the data and copy everything in a single folder (let ${DataFolder} be this folder)
2. move to the experiments folder:
   ```bash
   cd paper/experiments
   ```
3. run the following script (${DesiredFolder} is the folder where you want to store all the data used and produced during these experiments)
   ```bash
   bash moveData.sh ${DataFolder} ${DesiredFolder}
   ```

### Prerequisites
Before running the experiments, install the following prerequisites:
```bash
sudo apt-get update && sudo apt-get -y \
install build-essential g++ wget git cmake unzip zlib1g-dev \
python python3 python3-pip python3-biopython python3-biopython-sql \
python3-pysam python3-gffutils python-h5py python-scipy python-pysam \
emboss samtools cutadapt
```

### One-step script
Move to asgal repository root and run the following commands:
```bash
cd paper/experiments
bash runExperiments.sh ${DesiredFolder}
```

This script will run the full experiments, downloading and storing everything (**~60GB of space are required**) inside the specified folder:
- the folder _Tools_ will contain the tools used through the experiments ([Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml), [STAR](https://github.com/alexdobin/STAR), [SplAdder](http://raetschlab.org/suppl/spladder), and [AStalavista](http://sammeth.net/confluence/display/ASTA/Home))
- the folder _SimData_ will contain the simulated data
- the folder _AlComp_ will contain the output of _ASGAL_ aligner, _Hisat2_, and _STAR_ on simulated data
- the folder _SimEvents_ will contain the output of _ASGAL_ and _SplAdder_ on simulated data
- the folder _RealData_ will contain the output of _ASGAL_ and _SplAdder_ on real data
- **the folder _Results_ will contain a set of _CSV_ files summarizing the obtained results**

The one-step script is fully sequential and, on our server (Four 8-core Intel Xeon 2.30GHz processors and 256GB of RAM) it required ~10 days. As shown [below](#parallel), the experiments can be run in parallel.

##### Detailed explanation
If you want to run the experiments step by step, move to the experiments folder and
1. download _Hisat2_, _STAR_, _SplAdder_, and _AStalavista_
   ```bash
   bash downloadTools.sh ${DesiredFolder}
   ```
2. download simulated data (if you did not download them manually)
   ```bash
   bash downloadSimData.sh ${DesiredFolder}
   ```
3. set up the simulated data
   ```bash
   bash setupSimData.sh ${DesiredFolder}
   ```
3. run _ASGAL_ aligner, _Hisat2_, and _STAR_ on simulated data (_5M_ dataset)
   ```bash
   bash runAlignersComparison.sh ${DesiredFolder}
   ```
4. run _ASGAL_ and _SplAdder_ on simulated data (_5M_ and _10M_ datasets): this script will build the truth, create the reduced annotations, run _ASGAL_ full pipeline, _Hisat2_ + _SplAdder_ and _STAR_ + _SplAdder_ 
   ```bash
   bash runSimEvents.sh ${DesiredFolder}
   ```
5. download the real data (if you did not download them manually)
   ```bash
   bash downloadRealData.sh ${DesiredFolder}
   ```
6. run _ASGAL_ on real data: this script will build the truth, create the reduced annotations and run _ASGAL_ full pipeline
   ```bash
   bash runRealAsgal.sh ${DesiredFolder}
   ```
7. run _SplAdder_ on real data (4 different setups): this script will run _SplAdder_ using the reduced annotation **computed in the previous step** and the spliced alignments computed by _STAR_ and _Hisat2_
   ```bash
   bash runRealSplAdder.sh ${DesiredFolder}
   ```

###### <a name="parallel"></a>Parallel execution
The seven steps can be run in parallel reducing the time needed to complete the experiments (**the parallel execution has not been tested**):
1. run step 1
2. run step 2 and 6
3. run step 3
3. run step 4, 5, and 7
4. run step 8

##### Printing the results
To analyze the _CSV_ files and print the results, run the following commands:
- to print the accuracy and performance results of _ASGAL_ aligner, _Hisat2_, and _STAR_ run:
  ```bash
  bash ./AlComp/printResults.sh ${DesiredFolder}/Results/AlComp
  ```
- to print the accuracy and performance results of _ASGAL_ and _SplAdder_ on simulated data, run:
  ```bash
  bash ./SimEvents/printResults.sh ${DesiredFolder}/Results/SimEvents
  ```
- to print the accuracy and performance results of _ASGAL_ on real data, run:
  ```bash
  bash ./RealEvents/printAsgalResults.sh ${DesiredFolder}/Results/RealEvents/asgal
  ```
- to print the accuracy and performance results of _SplAdder_ on real data, run:
  ```bash
  bash ./RealEvents/printSpladderResults.sh ${DesiredFolder}/Results/RealEvents/spladder
  ```

### Printing our results
If you do not want to run the full experiments, we make available the _CSV_ files we computed and used in our experiments.
To print directly these results, run the following commands:
```bash
# move to asgal repository root
cd paper/experiments
tar xfz Results.tar.gz
bash ./AlComp/printResults.sh ./Results/AlComp
bash ./SimEvents/printResults.sh ./Results/SimEvents
bash ./RealEvents/printAsgalResults.sh ./Results/RealEvents/asgal
bash ./RealEvents/printSpladderResults.sh ./Results/RealEvents/spladder
```