# Use Ubtunu 20.04, using STAR to align.
# Adapted from the Cebola Lab's RNA-seq pipeline posted here:
# https://github.com/CebolaLab/RNA-seq
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /bin

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
        fastqc

# Install MultiQC
RUN pip install multiqc
