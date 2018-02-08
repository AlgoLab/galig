[Documentation](documentation) [Experiments](experiments)

**ASGAL** (**A**lternative **S**plicing **G**raph **AL**igner) is a
tool for detecting the alternative splicing events expressed in a
RNA-Seq sample with respect to a gene annotation. The 
_ASGAL_ approach consists of some  steps:
1. **splicing graph construction**: from the gene annotation, _ASGAL_ builds the splicing graph representing the gene structure that is implied by the set of input transcripts.
2. **splice-aware alignment**: _ASGAL_ aligns the RNA-Seq reads against the splicing
graph of the input gene. This procedure is tailored for such kind of alignments.
3. **detection of the alternative splicing events**: the spliced alignments are analyzed to detect the alternative splicing events that are induced by the reads in the sample. Moreover _ASGAL_ can report all events found, or only those that are not in the input annotation.

### Alternative Splicing Events
Actually, _ASGAL_ fully supports the following alternative splicing
events:
* _exon skipping_
* _alternative acceptor site_
* _alternative donor site_
* _intron retention_ (caused by the insertion of a new intron inside an
exon)

### Citation

If you use _ASGAL_, please cite its use as:

Luca Denti, Raffaella Rizzi, Stefano Beretta, Gianluca Della Vedova, Marco Previtali and Paola Bonizzoni.
_ASGAL: Aligning RNA-Seq Data to a Splicing Graph to Detect Novel Alternative Splicing Events_ (preprint available at [bioRxiv](https://doi.org/10.1101/260372))

### Install

_ASGAL_ is available <!--at conda-forge and--> at Docker Hub. A detailed installation walkthrough and the documentation is [available](documentation). 

### Example

The `example` directory contains a small dataset on the gene
[CG13375](http://www.ensembl.org/Drosophila_melanogaster/Gene/Summary?db=core;g=FBgn0040370;r=X:283186-294962)
of Drosophila Melanogaster.

If you already have [docker](https://www.docker.com) installed, you
can run _ASGAL_ on that sample with the commands
```
wget https://github.com/AlgoLab/galig/raw/master/example/input.tar.gz
tar xfz input.tar.gz
docker run -v "$PWD"/input:/data algolab/asgal
```

The running times will be a few seconds. Then you will find the file
[events.events](https://github.com/AlgoLab/galig/raw/master/example/events.events)
in the `input` directory. An extended explanation of this example
can be found <a href="http://asgal.algolab.eu/documentation#example"
target="_blank">here</a>.

### Contacts
If you have any question or you have any problem using the tool,
please contact [Luca
Denti](https://algolab.eu/people/luca-denti/). Alternatively, you can
use the [issue tracker](https://github.com/AlgoLab/galig/issues).

<!--
Given a gene annotation, the splicing graph is a graph where each
vertex is an exon and two vertices are linked if they are consecutive
in at least one transcript.
-->
