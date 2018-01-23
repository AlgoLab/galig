FROM ubuntu:17.10
RUN apt-get update && apt-get install -y \
    cmake \
    g++ \
    git-core \
    graphviz \
    python3-biopython \
    python3-gffutils \
    python3-graphviz \
    python3

VOLUME ["/data"]

RUN git clone https://github.com/AlgoLab/galig.git && \
    cd galig && \
    make prerequisites && \
    make

ENTRYPOINT ["/galig/asgal"]

CMD ["-g /data/DrosMel.BDGP6.chrX.fa", "-a /data/CG13375.gtf", "-s /data/CG13375.fasta", "-o /data/CG13375.events", "--allevents"]
