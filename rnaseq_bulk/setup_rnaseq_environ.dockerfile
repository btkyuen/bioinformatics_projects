# Use Ubtunu 20.04, using STAR to align.
# Adapted from the Cebola Lab's RNA-seq pipeline posted here:
# https://github.com/CebolaLab/RNA-seq
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN mkdir -p /home/apps

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
        samtools && \
    apt-get clean

# Install FastP
RUN cd /home/apps && \
    wget http://opengene.org/fastp/fastp && \
    chmod a+x ./fastp && \
    mv fastp /bin

# Install MultiQC
RUN pip install multiqc

# Install STAR aligner
# The STAR manual can be found here: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
RUN cd /home/apps/ && \
    wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz && \
    tar -xzf 2.7.11b.tar.gz && \
    rm 2.7.11b.tar.gz
RUN cp /home/apps/STAR-2.7.11b/bin/Linux_x86_64/STAR /bin && \
    rm -r /home/apps/STAR-2.7.11b

#RUN cd /home/apps && \
#    wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz && \
#    wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz && \
#    gunzip Homo_sapiens.GRCh38.113.gtf.gz Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Running genome generation step with less RAM, change this depending on the system configuration
#RUN STAR --runMode genomeGenerate \
#        --runThreadN 8 \
#        --genomeDir /home/apps/STAR/index \
#        --genomeFastaFiles /home/apps/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#        --sjdbGTFfile /home/apps/Homo_sapiens.GRCh38.113.gtf \
#        --sjdbOverhang 100 \
#        --limitGenomeGenerateRAM 14000000000 \
#        --genomeChrBinNbits 16 \
#        --genomeSAsparseD 2 

# Set environment variables for STAR
#RUN echo 'STAR_GENOME=/home/apps/STAR/index' >> /etc/environment && \
#    source /etc/environment


# Install deeptools
RUN pip3 install deeptools

# Install miniconda, then install salmon
RUN cd /home/apps && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /home/apps/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/home/apps/miniconda/bin":$PATH
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y -c bioconda salmon
RUN echo 'export PATH="/home/apps/miniconda/bin":$PATH' >> ~/.bashrc && \
    source ~/.bashrc

# Install GFF utilities for processing genome index for Salmon
RUN cd /home/apps && \
    wget http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.7.Linux_x86_64.tar.gz && \
    tar -xvf gffread-0.12.7.Linux_x86_64.tar.gz && \
    mv gffread-0.12.7.Linux_x86_64/gffread /bin && \
    rm -r gffread-0.12.7 && \
    rm gffread-0.12.7.tar.gz

# Generate Salmon index (moved to separate shell script)
#RUN mkdir -p /home/apps/salmon/index && \
#    gffread -w /home/apps/salmon/index/Homo_sapiens.GRCh38.dna.primary_assembly.transcripts.fa \
#            -g /home/apps/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#            /home/apps/Homo_sapiens.GRCh38.113.gtf && \
#    mv /home/apps/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai /home/apps/salmon/index/

# Set environment variables for Salmon
#RUN echo 'SALMON_GENOME=/home/apps/salmon/index/Homo_sapiens.GRCh38.dna.primary_assembly.transcripts.fa' >> /etc/environment && \
#    source /etc/environment

# Cleanup genome generation files
#RUN rm /home/apps/Homo_sapiens.GRCh38.113.gtf /home/apps/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Setup complete
RUN echo "Bulk RNAseq environment setup complete."