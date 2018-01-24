import sys, os

import gffutils

#Read gtf
def openGTF(gtfPath):
    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath),
                                 keep_order=True)
    except ValueError:
        gtf = gffutils.create_db(gtfPath,
                                 dbfn="{}.db".format(gtfPath),
                                 force=True, keep_order=True,
                                 disable_infer_genes=True,
                                 disable_infer_transcripts=True,
                                 merge_strategy='merge',
                                 sort_attribute_values=True)
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath), keep_order=True)
    return gtf

#Extract introns from gtf
def extractFromGTF(gtf):
    transcripts = []
    
    for g in gtf.features_of_type('gene'):
        for tr in gtf.children(g, featuretype='transcript', order_by='start'):
            transcript = []
            for ex in gtf.children(tr, featuretype='exon', order_by='start'):
                transcript += [ex.start]
                transcript += [ex.end]
            transcripts.append(transcript)
    return transcripts

def pairwise(iterable):
    a = iter(iterable)
    return zip(a, a)

def rewriteAnnotation(gtf, toRemove, outGTFPath):
    if all(toRemove):
        return
    if not(any(toRemove)):
        return
    outGTF = open(outGTFPath, 'w')
    for gene in gtf.features_of_type('gene'):
        outGTF.write(str(gene) + "\n")
        i = 0
        for tr in gtf.children(gene, featuretype='transcript', order_by='start'):
            if toRemove[i]:
                pass
            else:
                outGTF.write(str(tr) + "\n")
                for ex in gtf.children(tr, featuretype='exon', order_by='start'):
                    outGTF.write(str(ex) + "\n")
            i+=1
    outGTF.close()
    return

def main():
    gtfPath = sys.argv[1]
    eventsPath = sys.argv[2]
    outFold = sys.argv[3]

    gene = os.path.basename(gtfPath)[:-4]
    chrom = os.path.basename(os.path.dirname(gtfPath))

    if not os.path.isdir(os.path.join(outFold, chrom, gene)):
        os.makedirs(os.path.join(outFold, chrom, gene))
    gtf = openGTF(gtfPath)
    transcripts = extractFromGTF(gtf)

    counters = {'ES':0, 'A3':0, 'A5':0, 'IR':0}
    for line in open(eventsPath).readlines():
        info = line.strip("\n").split(" ")
        ev = info[1]

        uniqueName = ""
        if ev == 'ES' or ev == 'IR':
            counters[ev]+=1
            toRemove = [False]*len(transcripts)
            x,y = int(info[2]), int(info[3])
            i = 0
            for transcript in transcripts:
                for p1,p2 in pairwise(transcript[1:-1]):
                    if p1==min(x,y) and p2==max(x,y):
                        toRemove[i] = True
                i+=1
            uniqueName = ev + str(counters[ev]) + ".gtf"
            outGTFPath = os.path.join(outFold, chrom, gene, uniqueName)
            rewriteAnnotation(gtf, toRemove, outGTFPath)

        if ev == 'A3' or ev == 'A5':
            x1, x2, y = int(info[2]), int(info[3]), int(info[4])
            
            counters[ev]+=1
            toRemove = [False]*len(transcripts)
            i = 0
            for transcript in transcripts:
                for p1,p2 in pairwise(transcript[1:-1]):
                    if p1==min(x1,y) and p2==max(x1,y):
                        toRemove[i] = True
                i+=1
            if abs(y-x1)<abs(y-x2):
                uniqueName = ev + str(counters[ev]) + "e" + ".gtf"
            else:
                uniqueName = ev + str(counters[ev]) + "i" + ".gtf"
            outGTFPath = os.path.join(outFold, chrom, gene, uniqueName)
            rewriteAnnotation(gtf, toRemove, outGTFPath)
            x2 = int(info[3])
            toRemove = [False]*len(transcripts)
            i = 0
            for transcript in transcripts:
                for p1,p2 in pairwise(transcript[1:-1]):
                    if p1==min(x2,y) and p2==max(x2,y):
                        toRemove[i] = True
                i+=1
            if abs(y-x1)<abs(y-x2):
                uniqueName = ev + str(counters[ev]) + "i" + ".gtf"
            else:
                uniqueName = ev + str(counters[ev]) + "e" + ".gtf"
            outGTFPath = os.path.join(outFold, chrom, gene, uniqueName)
            rewriteAnnotation(gtf, toRemove, outGTFPath)

if __name__ == '__main__':
    main()
