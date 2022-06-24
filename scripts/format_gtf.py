import sys
import gffutils


def main():
    gtf_path = sys.argv[1]
    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtf_path), keep_order=True)
    except ValueError:
        gtf = gffutils.create_db(
            gtf_path,
            dbfn="{}.db".format(gtf_path),
            force=True,
            keep_order=True,
            disable_infer_genes=True,
            disable_infer_transcripts=False,
            merge_strategy="merge",
            sort_attribute_values=True,
        )

    for gene in gtf.features_of_type("gene"):
        if (
            len(list(gtf.children(gene, featuretype="transcript", order_by="start")))
            == 0
        ):
            continue
        print(gene)
        for transcript in gtf.children(
            gene, featuretype="transcript", order_by="start"
        ):
            print(transcript)
            for exon in gtf.children(transcript, featuretype="exon", order_by="start"):
                print(exon)


if __name__ == "__main__":
    main()
