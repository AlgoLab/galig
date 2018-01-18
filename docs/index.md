[Documentation](documentation) [Examples](examples) [Experiments](experiments)

**ASGAL** (**A**lternative **S**plicing **G**raph **AL**igner) is a
tool for detecting the alternative splicing events expressed in a
RNA-Seq sample with respect to a gene annotation. The **main idea**
behind ASGAL is the following one: the alternative splicing events can
be detected by aligning the RNA-Seq reads against the splicing graph
of the gene.

ASGAL approach can be divided in two steps:
1. **splice-aware alignment** of the RNA-Seq reads against the splicing
graph of the input gene
2. **detection of the alternative splicing events** expressed in the
sample with respect to the input annotation

The main difference between ASGAL and the other tools for the
detection of alternative splicing events is that ASGAL considers each
annotated transcript as an expressed isoform and detects all the
differences between them and the input RNA-Seq reads. For this reason,
we say that:
> ASGAL detects the alternative splicing events expressed
> in a RNA-Seq sample with respect to a given gene annotation

### Alternative Splicing Events
Actually, ASGAL fully supports the following alternative splicing
events:
* _exon skipping_
* _alternative acceptor site_
* _alternative donor site_
* _intron retention_ (caused by the insertion of a new intron inside an
exon)

We are working on extending our approach in order to support:
* _(small) cassette exon_
* (more complex but less probable) events arising from the combination
  of two or more events

ASGAL will never be able to support, since it aligns the reads to the
splicing graph and it cannot align to the introns of the gene:
* _(long) cassette exon_
* _intron retention_ (caused by the union of two exons)

--- figure ---

<!--
### Citations
-->

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
