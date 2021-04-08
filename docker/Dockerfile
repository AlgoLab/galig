FROM ubuntu:20.04

## for apt to be noninteractive
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt update -qqy
RUN apt-get install -qqy apt-utils

RUN echo 'tzdata tzdata/Areas select Etc' | debconf-set-selections; \
    echo 'tzdata tzdata/Zones/Etc select UTC' | debconf-set-selections; \
    apt-get install -qqy --no-install-recommends tzdata

## preesed tzdata, update package index, upgrade packages and install needed software
RUN apt-get install -qqy \
    build-essential \
    ca-certificates \
    cmake \
    dirmngr \
    git \
    make \
    python3 \
    python3-biopython \
    python3-biopython-sql \
    python3-pandas \
    python3-pip \
    python3-pysam \
    python3-setuptools \
    samtools \
    util-linux \
    wget \
    zlib1g-dev
RUN apt-get clean
RUN pip3 install gffutils

RUN git clone --depth=1 --recursive https://github.com/AlgoLab/galig.git
RUN cd galig ; make prerequisites
RUN cd galig ; make

VOLUME ["/data"]
COPY asgal-docker.sh /galig/asgal-docker.sh
ENV PATH=$PATH:/galig/asgal
ENTRYPOINT ["/galig/asgal-docker.sh"]
