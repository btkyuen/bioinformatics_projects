# Use Ubtunu 20.04, using STAR to align.
# Adapted from the NCI's DNA Sequence Variant Calling Pipeline
# found here: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Set up Ubuntu requirements for software
# Install FastQC
RUN apt-get update && \
    apt-get install -y -q --no-install-recommends\
        wget \
        gfortran \
        build-essential \
        software-properties-common \
        libcurl4-openssl-dev \
        libssl-dev \
        libfontconfig1-dev \
        libxml2-dev \
        python3.9 \
        python3-pip \
        python3-setuptools \
        python3-dev \
        fastqc \ 
        git \
        openjdk-17-jdk

# Install MultiQC
RUN pip install multiqc

# Install Burrows-Wheeler Aligner (bwa)
# Install both bwa-mem and bwa-aln, to handle a range of sequences
RUN mkdir -p ~/apps && \
    cd ~/apps && \
    git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make 
RUN cp ~/apps/bwa/bwa /bin

# Install Picard tools, and add a shortcut to the environment
RUN mkdir -p ~/apps/picard && \
    cd ~/apps/picard && \
    wget https://github.com/broadinstitute/picard/releases/download/3.3.0/picard.jar
RUN echo PICARD=/root/apps/picard/picard.jar >> /etc/environment && \
    source /etc/environment

