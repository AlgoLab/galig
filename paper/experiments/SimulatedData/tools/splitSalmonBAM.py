#!/usr/bin/python3

import argparse
import os
import sys
import copy
import random

import pysam
import sys
import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip

def build_gene_dict(annotation_filename):
    tr_gene_dict = {}
    db_filename = "{}.db".format(annotation_filename)
    if not os.path.exists(db_filename):
        gffutils.create_db(annotation_filename, dbfn=db_filename)
    gtf_file = gffutils.FeatureDB(db_filename, keep_order = True)
    for gene in gtf_file.features_of_type('gene'):
        for tr in gtf_file.children(gene, featuretype = "transcript"):
            tr_gene_dict[tr.id] = gene.id
    return tr_gene_dict


def split_bam(aln_filename, gene_chr_dict, gene_list, out_dir, unmapped_reads):
    out_files = {}
    input_mode = 'rb' if aln_filename.endswith('bam') else 'r'
    with pysam.AlignmentFile(aln_filename, input_mode) as aln_file:
        for aln in aln_file:        
            tr_name = aln_file.get_reference_name(aln.reference_id)
            if tr_name not in gene_chr_dict:
                continue
            if gene_chr_dict[tr_name] in gene_list:
                outf = gene_chr_dict[tr_name]
                if gene_chr_dict[tr_name] not in out_files:
                    output_filename = "{}/{}.fa.gz".format(out_dir, outf)
                    file_handle = gzip.open(output_filename, "wt")
                    out_files[outf] = file_handle
                fasta_seq = SeqRecord(Seq(aln.query_sequence), id=aln.query_name, name=aln.query_name, description="From Salmon alignment file")
                SeqIO.write(fasta_seq, out_files[outf], "fasta")
                if fasta_seq.id in unmapped_reads:
                    SeqIO.write(unmapped_reads[fasta_seq.id], out_files[outf], "fasta")

def parse_unmapped_file(unmapped_filename, sample_filenames):
    unmapped_tags = set()
    with open(unmapped_filename) as unmapped_file:
        for line in unmapped_file:
            tag, unmap_type = line.strip('\n').split(' ')
            if unmap_type in ['m1', 'm2']:
                file_suff = "2" if unmap_type == 'm1' else "1"
                unmapped_tags.add((tag, file_suff))
    print("Parsing samples to retrieve reads")
    unmapped_reads = {}
    for sample_fn in sample_filenames:
        print("Parsing {}".format(sample_fn))
        file_name = os.path.basename(sample_fn)
        file_suff = file_name[file_name.find('_')+1:-len('.fastq.gz')]
        file_handle = gzip.open(sample_fn, "rt")
        for record in SeqIO.parse(file_handle, "fastq"):
            test_elem = (record.id, file_suff)
            if test_elem in unmapped_tags:
                record_copy = copy.deepcopy(record)
                record_copy.id="{}/{}".format(record_copy.id, file_suff)
                unmapped_reads[record.id] = record_copy
    return unmapped_reads

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('alignment_filename', metavar='BAM/SAM FILE', help='Salmon output BAM/SAM')
    parser.add_argument('annotation_filename', metavar='GTF FILE', help='Annotation file')
    parser.add_argument('unmapped_filename', metavar='UNMAPPED FILE', help='Unmapped file from Salmon')
    parser.add_argument('sample_filenames', metavar='SAMPLE FILES', nargs=2, help='Reads sample')
    parser.add_argument('-g', '--gene-list', metavar='GENE FILE', dest='genes',
                        help='File with a list of target genes', required = 'true')
    parser.add_argument('-o', '--output-dir', dest='output_dir', help='Output Directory',
                        default='.')
    args = parser.parse_args()

    out_dir = os.path.abspath(args.output_dir)

    if not os.path.isfile(args.alignment_filename):
        print("Error: {} is not a file.".format(args.alignment_filename))
        print("Aborting...")
        sys.exit(1)
    
    if os.path.exists(out_dir) and os.path.isfile(out_dir):
        print("Error: {} is a file.".format(out_dir))
        print("Aborting...")
        sys.exit(1)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    print("Parsing GTF annotation")
    gene_chr_dict = build_gene_dict(args.annotation_filename)

    print("Retrieving unmapped reads")
    unmapped_reads = parse_unmapped_file(args.unmapped_filename, args.sample_filenames)
    print("unmapped reads that will be remapped: {}".format(len(unmapped_reads)))

    print("Split BAM")
    gene_list = list()
    num_genes = 0
    tot_genes = 0

    with open(args.genes, "r") as gf:
        for g in gf.readlines():
            if num_genes == 930:
                print("Parsing BAM alignments")    
                split_bam(args.alignment_filename, gene_chr_dict, gene_list, out_dir, unmapped_reads)
                num_genes = 0
                gene_list = list()
                print("Genes completed: {}".format(tot_genes))
            else:
                num_genes += 1
                gene_list.append(g.rstrip())
            tot_genes += 1
        if num_genes > 0:
            print("Parsing BAM alignments")    
            split_bam(args.alignment_filename, gene_chr_dict, gene_list, out_dir, unmapped_reads)
    return

if __name__ == "__main__":
    main()
