import gffutils
from Bio import SeqIO

class Annotation:
    def __init__(self, genomic_path, gtf_path):
        genomic = SeqIO.index(genomic_path, "fasta")
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
        T = "|"
        edges = []
        added_exons = {}
        ex_index = 1

        self.gene = ""
        self.transcripts = {}
        self.exons = {}
        self.introns = {}
        self.exons_name = [""]
        self.exons_length = [-1]
        self.strand = True
        for gene in gtf.features_of_type('gene'):
            chromo = gene.seqid
            sequence = genomic[chromo]
            self.strand = gene.strand == '+'
            self.gene = gene.id
            for tr in gtf.children(gene, featuretype='transcript', order_by='start'):
                t_id = tr.id
                self.transcripts.update({t_id : tr.start})
                prev_start = -1
                Es = []
                Is = []
                last_id = -1
                for ex in gtf.children(tr, featuretype='exon', order_by='start'):
                    start = ex.start
                    end = ex.end
                    ex_posid = str(start) + "_" + str(end)
                    Es.append([start, end])
                    if prev_start != -1:
                        Is.append([prev_start, start-1])
                    prev_start = end+1
                    if ex_posid not in added_exons:
                        ex_seq = sequence[ex.start-1:ex.end-1+1].seq
                        T += ex_seq + "|"
                        self.exons_length.append(len(ex_seq))
                        self.exons_name.append(ex.attributes['exon_id'][0])
                        added_exons.update({ex_posid:ex_index})
                        ex_index+=1
                    if last_id != -1:
                        edges.append((last_id, ex_posid))
                    last_id = ex_posid
                self.exons.update({t_id : Es})
                self.introns.update({t_id : Is})

        self.transcripts_length = {}
        for tr_id,Es in self.exons.items():
            l = 0
            for (s,e) in Es:
                l+=e-s
            self.transcripts_length.update({tr_id:l})

        self.BV = BitVector(T)

        self.adj_matrix = []
        for i in range(0, ex_index+1):
            row = []
            for j in range(0, ex_index+1):
                row.append(0)
            self.adj_matrix.append(row)
        for edge in edges:
            self.adj_matrix[added_exons[edge[0]]][added_exons[edge[1]]] = 1

    def getExonLength(self, i):
        return self.exons_length[i]

    def getExonName(self, i):
        return self.exons_name[i]
    
    def contain(self, i, j):
        return self.adj_matrix[i][j]

    def rank(self, i):
        return self.BV.rank(i)

    def select(self, i):
        return self.BV.select(i)
    
class BitVector:
    '''
    Il bit vector indicizza da 0. 
    Presenta un 1 dove nel testo di trova un "|"
    '''
    def __init__(self, text):
        self.text = text
        
        self.length = len(text)
        
        self.bit_vector = []
        self.ones = []
        i = 1
        for c in self.text:
            if c == "|":
                self.bit_vector.append(1)
                self.ones.append(i)
            else:
                self.bit_vector.append(0)
            i += 1

    def rank(self, i):
        n = 0
        for elem in self.ones:
            if elem<=int(i):
                n+=1
            else:
                break
        return n
    
    def select(self, i):
        iterator = 0
        for elem in self.ones:
            iterator+=1
            if iterator == i:
                return elem

    def __str__(self):
        s = ""
        for elem in self.bit_vector:
            s += str(elem)
        return s

if __name__ == '__main__':
    pass
