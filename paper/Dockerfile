FROM ubuntu:18.04

RUN apt update && \
    apt install -y build-essential \
    cmake \
    make \
    wget \
    git \
    python3-pip \
    samtools \
    zlib1g-dev \
    python3 \
    python3-pip \
    python3-setuptools \
    python3-biopython \
    python3-biopython-sql \
    python3-pysam \
    python3-pandas \
    python3-setuptools ; \
    pip3 install gffutils

VOLUME ["/data"]

RUN git clone --recursive https://github.com/AlgoLab/galig.git ; \
    cd galig ; \
    make prerequisites ; \
    make

ENTRYPOINT ["/galig/asgal-docker.sh"]
