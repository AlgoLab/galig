import sys, os, re

from Bio import SeqIO
import gffutils

def main():
    gtf_path = sys.argv[1]
    fa_path = sys.argv[2]
    sam_path = sys.argv[3]

    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtf_path),
                                 keep_order=True)
    except ValueError:
        gtf = gffutils.create_db(gtf_path,
                                 dbfn="{}.db".format(gtf_path),
                                 force=True, keep_order=True,
                                 disable_infer_genes=True,
                                 disable_infer_transcripts=True,
                                 merge_strategy='merge',
                                 sort_attribute_values=True)

    transcripts = {}
    exons = {}
    introns = {}
    for gene in gtf.features_of_type('gene'):
        STRAND = gene.strand == '+'
        gene = gene.id
        for tr in gtf.children(gene, featuretype='transcript', order_by='start'):
            t_id = tr.id
            transcripts.update({t_id : [tr.start, tr.end]})
            prev_start = -1
            Es = []
            Is = []
            for ex in gtf.children(tr, featuretype='exon', order_by='start'):
                start = ex.start
                end = ex.end
                Es.append([start, end])
                if prev_start != -1:
                    Is.append([prev_start, start-1])
                prev_start = end+1
            exons.update({t_id : Es})
            introns.update({t_id : Is})

    TP = 0
    FP = 0
    FN = 0

    found_ids = []
    with open(sam_path, 'r') as sam:
        for line in sam:
            if line[0] != '@':
                spl_line = line.split('\t')
                identifier, strand, pos1, cigar, pos2 = spl_line[0], int(spl_line[1])==0, int(spl_line[3]), spl_line[5], int(spl_line[12].split(":")[2])
                found_ids.append(identifier)
                regex = re.search("(.+)\/(.+);mate1:(.+)-.+;mate2:(.+)-.+$$", identifier)
                read_n, transcript, start1, start2 = regex.group(1), regex.group(2), int(regex.group(3)), int(regex.group(4))

                if STRAND:
                    initial_clipping = 0
                    try:
                        initial_clipping = int(re.search("^[0-9]+S", cigar).group(0)[:-1])
                    except AttributeError:
                        initial_clipping = 0
                    realPos = 1-initial_clipping
                    for E in exons[transcript]:
                        if E[1] > pos1:
                            realPos+=pos1-E[0]
                            break
                        else:
                            realPos+=E[1]-E[0]+1
                    if strand:
                        if start1 == realPos:
                            TP+=1
                        else:
                            print("+", read_n, start1, realPos, cigar)
                            FP+=1
                    else:
                        if start2 == realPos:
                            TP+=1
                        else:
                            print("-", read_n, start2, realPos, cigar)
                            FP+=1
                else:
                    final_clipping = 0
                    try:
                        final_clipping = int(re.search("[0-9]+S$", cigar).group(0)[:-1])
                    except AttributeError:
                        final_clipping = 0
                    realPos = 1-final_clipping
                    for E in reversed(exons[transcript]):
                        if E[0] < pos2:
                            realPos+=E[1]-pos2
                            break
                        else:
                            realPos+=E[1]-E[0]+1
                    realPos += 7
                    #print(read_n, start1, start2, realPos, cigar)
                    if strand:
                        if start2 == realPos:
                            TP+=1
                        else:
                            print("+", read_n, start2, realPos, cigar)
                            FP+=1
                    else:
                        if start1 == realPos:
                            TP+=1
                        else:
                            print("-", read_n, start1, realPos, cigar)
                            FP+=1

    real = 0
    for record in SeqIO.parse(fa_path, "fasta"):
        real+=1
        if record.id not in found_ids:
            FN+=1

    P = TP/(TP+FP)
    R = TP/(TP+FN)
    print(gene, real, len(found_ids), TP, FP, FN, P, R)

if __name__ == '__main__':
    main()
