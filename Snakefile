configfile: "config.yaml"

import re
from os.path import join as pjoin

FA = config["fa"]
GTF = config["gtf"]
TR = config["tr"]

FQ1 = config["fq1"]
FQ2 = config["fq2"]

ODIR = config["odir"]

THREADS = config["threads"]

ASGAL_DIR = config["galig"]

genes = {}
for line in open(GTF):
    if line.startswith("#"):
        continue
    line = line.strip("\n").split("\t")
    chrom = line[0]
    if line[2] == "gene":
        gene = re.match("gene_id \"([A-Za-z0-9\.]+)\";", line[-1]).group(1)
        genes[gene] = chrom
print(genes)

rule run:
    input:
        expand(pjoin(ODIR, "{gene}", "ASGAL", "events.wpsi.csv"),
                gene = genes.keys())

rule salmon_index:
    input:
        fa = TR
    output:
        index = directory(pjoin(ODIR, "salmon-index"))
    threads: THREADS
    shell:
        """
        salmon index -p {threads} -t {input.fa} -i {output.index}
        """

rule salmon_quant:
    input:
        fq1 = FQ1,
        fq2 = FQ2,
        index = pjoin(ODIR, "salmon-index")
    output:
        bam = pjoin(ODIR, "salmon-out", "aligns.bam"),
        txt = pjoin(ODIR, "salmon-out", "aux_info", "unmapped_names.txt")
    params:
        odir = pjoin(ODIR, "salmon-out"),
        sam = pjoin(ODIR, "salmon-out", "aligns.sam")
    threads: THREADS
    shell:
        """
        salmon quant -p {threads} -i {input.index} -l A -1 {input.fq1} -2 {input.fq2} -o {params.odir} --validateMappings --writeMappings={params.sam} --writeUnmappedNames
        samtools view -bS {params.sam} | samtools sort > {output.bam}
        samtools index {output.bam}
        """

rule split_reads:
    input:
        gtf = GTF,
        bam = pjoin(ODIR, "salmon-out", "aligns.bam")
    output:
        fq = pjoin(ODIR, "{gene}", "reads.fq")
    params:
        bed = pjoin(ODIR, "{gene}", "transcripts.bed")
    threads: 1
    shell:
        """
        grep {wildcards.gene} {input.gtf} | grep -P "\\ttranscript\\t" | egrep -o 'transcript_id "[0-9A-Za-z\.]+";' | cut -f 2 -d'"' | cut -f1 -d'.' | sort -u | awk '{{ print $0"\\t"1"\\t"4294967295 }}'> {params.bed}
        samtools view -b {input.bam} -L {params.bed} | samtools fastq > {output.fq}
        """

rule get_unmapped:
    input:
        fq1 = FQ1,
        fq2 = FQ2,
        txt = pjoin(ODIR, "salmon-out", "aux_info", "unmapped_names.txt")
    output:
        fq = pjoin(ODIR, "salmon-out", "unmapped_names.fq")
    params:
        txt = pjoin(ODIR, "salmon-out", "aux_info", "unmapped_names.one-col.txt")
    threads: 1
    shell:
        """
        cut -f1 -d' ' {input.txt} > {params.txt}
        grep -A 1 -f {params.txt} {input.fq1} | seqtk seq -F '#' > {output.fq}
        grep -A 1 -f {params.txt} {input.fq2} | seqtk seq -F '#' >> {output.fq}
        """

rule complete_fq:
    input:
        fq = pjoin(ODIR, "{gene}", "reads.fq"),
        unmapped_fq = pjoin(ODIR, "salmon-out", "unmapped_names.fq")
    output:
        fq = pjoin(ODIR, "{gene}", "reads.wunmapped.fq")
    threads: 1
    shell:
        """
        cat {input.fq} {input.unmapped_fq} > {output.fq}
        """

rule split_reference:
    input:
        fa = FA
    output:
        fa = pjoin(ODIR, "chroms", "{chrom}.fa")
    threads: 1
    shell:
        """
        samtools faidx {input.fa} {wildcards.chrom} > {output.fa}
        """

rule split_annotation:
    input:
        gtf = GTF
    output:
        gtf = pjoin(ODIR, "{gene}", "annotation.gtf")
    threads: 1
    shell:
        """
        grep {wildcards.gene} {input.gtf} > {output.gtf}
        """

rule asgal:
    input:
        fa = lambda wildcards: pjoin(ODIR, "chroms", f"{genes[wildcards.gene]}.fa"),
        gtf = pjoin(ODIR, "{gene}", "annotation.gtf"),
        fq = pjoin(ODIR, "{gene}", "reads.wunmapped.fq")
    output:
        sam = pjoin(ODIR, "{gene}", "ASGAL", "aligns.sam"),
        csv = pjoin(ODIR, "{gene}", "ASGAL", "events.csv")
    params:
        mem = pjoin(ODIR, "{gene}", "ASGAL", "aligns.mem")
    threads: 1
    shell:
        """
        {ASGAL_DIR}/bin/SpliceAwareAligner -g {input.fa} -a {input.gtf} -s {input.fq} -o {params.mem}
        python3 {ASGAL_DIR}/scripts/formatSAM.py -m {params.mem} -g {input.fa} -a {input.gtf} -o {output.sam}
        python3 {ASGAL_DIR}/scripts/detectEvents.py --allevents -g {input.fa} -a {input.gtf} -m {params.mem} -o {output.csv}
        """

rule sam2bam:
    input:
        "{f}.sam"
    output:
        "{f}.bam"
    threads: 1
    shell:
        """
        samtools view -bS {input} | samtools sort > {output}
        samtools index {output}
        """

rule compute_psi:
    input:
        gtf = pjoin(ODIR, "{gene}", "annotation.gtf"),
        bam = pjoin(ODIR, "{gene}", "ASGAL", "aligns.bam"),
        csv = pjoin(ODIR, "{gene}", "ASGAL", "events.csv")
    output:
        csv = pjoin(ODIR, "{gene}", "ASGAL", "events.wpsi.csv")
    shell:
        """
        python3 {ASGAL_DIR}/scripts/compute_psi.py {input.gtf} {input.bam} {input.csv} > {output.csv}
        """
