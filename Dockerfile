# Base Image
FROM biocontainers/biocontainers:latest

# Metadata
LABEL base.image="biocontainers:latest"
LABEL version="2"
LABEL software="seqcheck-deps"
LABEL software.version="0.1"
LABEL description="dependencies for mini-pipeline to sanity check new sequencing data on Synapse"
LABEL website="https://github.com/jaeddy/seqcheck"
LABEL documentation="https://github.com/jaeddy/seqcheck"
LABEL license="https://github.com/jaeddy/seqcheck"
LABEL tags="Genomics"

# Maintainer
MAINTAINER James Eddy <james.a.eddy@gmail.com>

# set up packages
RUN conda config --add channels conda-forge && \
    conda config --add channels defaults && \
    conda config --add channels r && \
    conda config --add channels bioconda
RUN conda install fastqc=0.11.5 hisat2=2.0.5 salmon=0.8.2 multiqc=0.9.1

USER root

# set version here to minimize need for edits below
# ENV VERSION=v1.6.2
ENV BRANCH=develop
ENV VERSION=6ba6a3ebde81fe8ed4d0c231ab42c613aa03334f

ENV PACKAGES python-dev git python-setuptools python-pip

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES}

RUN git clone -b ${BRANCH} git://github.com/Sage-Bionetworks/synapsePythonClient.git && \
    cd synapsePythonClient && \
    git checkout ${VERSION} && \
    python setup.py develop
