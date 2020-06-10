FROM ubuntu:18.04

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
    wget \
    zlib1g-dev
RUN apt-get clean
RUN pip3 install gffutils


# Install gosu
ENV GOSU_VERSION 1.12
RUN set -eux; \
# save list of currently installed packages for later so we can clean up
	dpkgArch="$(dpkg --print-architecture | awk -F- '{ print $NF }')"; \
	wget -O /usr/local/bin/gosu "https://github.com/tianon/gosu/releases/download/$GOSU_VERSION/gosu-$dpkgArch"; \
	wget -O /usr/local/bin/gosu.asc "https://github.com/tianon/gosu/releases/download/$GOSU_VERSION/gosu-$dpkgArch.asc"; \
	export GNUPGHOME="$(mktemp -d)"; \
	gpg --batch --keyserver hkps://keys.openpgp.org --recv-keys B42F6819007F00F88E364FD4036A9C25BF357DD4; \
	gpg --batch --verify /usr/local/bin/gosu.asc /usr/local/bin/gosu; \
	command -v gpgconf && gpgconf --kill all || :; \
	rm -rf "$GNUPGHOME" /usr/local/bin/gosu.asc; \
	\
# clean up fetch dependencies
	chmod +x /usr/local/bin/gosu; \
# verify that the binary works
	gosu --version; \
	gosu nobody true

RUN git clone --recursive https://github.com/AlgoLab/galig.git
RUN cd galig ; make prerequisites
RUN cd galig ; make

VOLUME ["/data"]
COPY asgal-docker.sh /galig/asgal-docker.sh
ENTRYPOINT ["/galig/asgal-docker.sh"]
