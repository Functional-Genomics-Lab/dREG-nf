# syntax=docker/dockerfile:1
FROM ubuntu:22.04

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    r-base-core \
    # libssl1.0.0 \
    libssl-dev \
    bzip2 \
    make \
    libcurl4-gnutls-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ADD --checksum=sha256:cbe2f6911e1455e4386160106f796eb54c6e2f8456467d4557f93ee135c66ac8 \
    https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2 \
    bedops_linux.tar.bz2

RUN tar jxvf bedops_linux.tar.bz2
RUN ln -s /usr/share/dotnet/dotnet /usr/bin/dotnet

ADD --checksum=sha256:ea508d763dc27d30d1a34f6a38d8da31ea67797d7b4971e8c10677bc48539565 \
    https://github.com/Danko-Lab/dREG.git \
    dREG

WORKDIR dREG

RUN make R_dependencies
RUN make dreg

ENV PATH dREG/:$PATH:/usr/local/bin
