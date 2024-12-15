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

# Install FastP
RUN wget http://opengene.org/fastp/fastp && \
    chmod a+x ./fastp

# Install MultiQC
RUN pip install multiqc

# Install STAR aligner
# The STAR manual can be found here: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
RUN wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz && \
    tar -xzf 2.7.11b.tar.gz && \
    rm 2.7.11b.tar.gz
RUN cp /bin/Linux_x86_64/STAR /bin/

# Get human genome files and start building indices
# Because read length could vary, setting the sjdbOverhang setting to 100bp
RUN mkdir -p /home/apps/STAR/index && cd /home/apps/STAR/index
RUN wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
RUN wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
RUN gunzip Homo_sapiens.GRCh38.113.gtf.gz Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# Running genome generation step with less RAM, change this depending on the system configuration
RUN STAR --runMode genomeGenerate --runThreadN 8 --genomeDir /home/apps/STAR/index --genomeFastaFiles /home/apps/STAR/index/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --sjdbGTFfile /home/apps/STAR/index/Homo_sapiens.GRCh38.113.gtf --sjdbOverhang 100 --limitGenomeGenerateRAM 10000000000





