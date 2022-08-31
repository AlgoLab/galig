import sys
import gffutils


def open_gtf(gtf_path):
    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtf_path), keep_order=True)
    except ValueError:
        gtf = gffutils.create_db(
            gtf_path,
            dbfn="{}.db".format(gtf_path),
            force=True,
            keep_order=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True,
            merge_strategy="merge",
            sort_attribute_values=True,
        )
    return gtf


def main():
    gtf_path = sys.argv[1]
    gtf = open_gtf(gtf_path)

    for gene in gtf.features_of_type("gene"):
        print(gene)
        fake_exons = set()
        for transcript in gtf.children(
            gene, featuretype="transcript", order_by="start"
        ):
            print(transcript)
            exons = [
                exon
                for exon in gtf.children(
                    transcript, featuretype="exon", order_by="start"
                )
            ]
            for exon in exons:
                print(exon)
            for i in range(len(exons) - 1):
                l = exons[i + 1].start - 1 - (exons[i].end + 1) + 1
                if l < 100 and l % 3 == 0:
                    fake_exons.add((exons[i].start, exons[i + 1].end))
        i = 1
        for (s, e) in fake_exons:
            print(
                gene.chrom,
                gene.source,
                "transcript",
                s,
                e,
                ".",
                gene.strand,
                ".",
                f'gene_id "{gene.attributes["gene_id"][0]}"; transcript_id "{gene.attributes["gene_id"][0]}_fake_transcript_{i}";',
                sep="\t",
            )
            print(
                gene.chrom,
                gene.source,
                "exon",
                s,
                e,
                ".",
                gene.strand,
                ".",
                f'gene_id "{gene.attributes["gene_id"][0]}"; transcript_id "{gene.attributes["gene_id"][0]}_fake_transcript_{i}"; exon_id "fake_exon_{i}"',
                sep="\t",
            )
            i += 1


if __name__ == "__main__":
    main()
