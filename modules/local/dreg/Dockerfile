# syntax=docker/dockerfile:1
# Ideal would be 3.6.0 like they used for Bidirectional-Flow
# https://github.com/Dowell-Lab/Bidirectional-Flow/blob/0e83e678b3321d9986e7eba9796d9c4828c01147/README.md?plain=1#L71
FROM rocker/tidyverse:3.6.3-ubuntu18.04
LABEL org.opencontainers.image.source https://github.com/Functional-Genomics-Lab/dREG-nf

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    # libssl1.0.0 \
    libssl-dev \
    bzip2 \
    make \
    libcurl4-gnutls-dev \
    libpcre3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ADD --checksum=sha256:cbe2f6911e1455e4386160106f796eb54c6e2f8456467d4557f93ee135c66ac8 \
    https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2 \
    bedops_linux.tar.bz2

ADD --checksum=sha256:2fc5fdd4dcf26ddf44427318925a69a743e97cb88ccc4964c314fcf389403584 \
    https://cran.r-project.org/src/contrib/Archive/rphast/rphast_1.6.9.tar.gz \
    rphast.tar.gz

RUN R CMD INSTALL rphast.tar.gz

ADD --checksum=sha256:f4b88920419eb0a89d0bc5744af0416d92d112988702dc726882394128a8754d \
    https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-14.tar.gz \
    randomForest.tar.gz

RUN R CMD INSTALL randomForest.tar.gz

ADD https://github.com/Danko-Lab/dREG.git dREG

WORKDIR dREG

RUN install2.r --error --skipinstalled --ncpus -1 \
    --repos http://cran.us.r-project.org --repos getOption \
    data.table \
    e1071 \
    mvtnorm \
    parallel \
    rmutil \
    snowfall \
    && rm -rf /tmp/downloaded_packages \
    && strip /usr/local/lib/R/site-library/*/libs/*.so

RUN installGithub.r andrelmartins/bigWig/bigWig

RUN R CMD INSTALL dREG --clean

ENV PATH dREG/:$PATH:/usr/local/bin

# RUN tar jxvf bedops_linux.tar.bz2 -C ./bedops
# RUN cp -r bedops/bin/ /usr/bin/bedops
