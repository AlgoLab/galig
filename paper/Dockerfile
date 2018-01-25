FROM ubuntu:17.10
RUN apt-get update && apt-get install -y \
    cmake \
    g++ \
    git-core \
    graphviz \
    python3-biopython \
    python3-biopython-sql \
    python3-gffutils \
    python3-graphviz \
    python3

VOLUME ["/data"]

RUN git clone https://github.com/AlgoLab/galig.git && \
    cd galig && \
    make prerequisites && \
    make

ENTRYPOINT ["/galig/asgal"]

CMD ["-g /data/genome.fa", "-a /data/annotation.gtf", "-s /data/reads.fasta", "-o /data/events"]
