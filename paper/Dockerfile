FROM ubuntu:18.04
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    wget \
    unzip \
    git \
    python3 \
    python3-pip \
    python3-biopython \
    python3-biopython-sql \
    python3-pysam \
    zlib1g-dev \
    libboost1.65-all-dev \
    curl \
    liblzma-dev \
    libjemalloc-dev \
    libjemalloc1 \
    libghc-bzlib-dev \
    samtools ; \
    pip3 install gffutils

VOLUME ["/data"]

RUN git clone --recursive https://github.com/AlgoLab/galig.git ; \
    cd galig ; \
    make prerequisites ; \
    make

ENTRYPOINT ["/galig/asgal"]

CMD ["-g /data/genome.fa", "-a /data/annotation.gtf", "-s /data/reads.fasta", "-o /data/events"]
