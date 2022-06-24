import sys
import pysam

def main():
    gtf_path = sys.argv[1]
    bam_path = sys.argv[2]
    csv_path = sys.argv[3]

    chrom = open(gtf_path).readline().split("\t")[0]

    bam = pysam.AlignmentFile(bam_path, "rb")

    for line in open(csv_path):
        if line.startswith("Type"):
            print("Type,Start,End,ExclusionSupport,InclusionSupport,PSI,Transcripts")
            continue
        etype, s, e, exclusion_w, Ts = line.strip("\n").split(",")
        s,e = int(s), int(e)
        exclusion_w = int(exclusion_w)
        all_w = 0
        for al in bam.fetch(chrom, s, e):
            if etype == "IR":
                all_w += 1
            elif al.get_cigar_stats()[1][3] > 0:
                # we have Ns
                all_w += 1
        inclusion_w = all_w - exclusion_w

        psi = round(inclusion_w/(inclusion_w+exclusion_w), 2) if inclusion_w+exclusion_w != 0 else 0
        print(etype, s, e, exclusion_w, inclusion_w, psi, Ts, sep=",")


if __name__ == "__main__":
    main()
