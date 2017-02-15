# Base Image
FROM biocontainers/biocontainers:latest

# Metadata
LABEL base.image="biocontainers:latest"
LABEL version="1"
LABEL software="seqcheck-deps"
LABEL software.version="0.1"
LABEL description="dependencies for mini-pipeline to sanity check new sequencing data on Synapse"
LABEL website="https://github.com/jaeddy/seqcheck"
LABEL documentation="https://github.com/jaeddy/seqcheck"
LABEL license="https://github.com/jaeddy/seqcheck"
LABEL tags="Genomics"

# Maintainer
MAINTAINER James Eddy <james.a.eddy@gmail.com>

RUN conda install synapseclient=1.5 fastqc=0.11.5 hisat2=2.0.5 salmon=0.7.2 
