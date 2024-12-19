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
        fastqc

# Install FastP
RUN cd /home/apps && \
    wget http://opengene.org/fastp/fastp && \
    chmod a+x ./fastp && \
    mv fastp /bin

# Install MultiQC
RUN pip install multiqc

# Install STAR aligner
# The STAR manual can be found here: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
RUN mkdir -p /home/apps/STAR/index && \
    cd /home/apps/ && \
    wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz && \
    tar -xzf 2.7.11b.tar.gz && \
    rm 2.7.11b.tar.gz
RUN cp /home/apps/STAR-2.7.11b/bin/Linux_x86_64/STAR /bin && \
    rm -r /home/apps/STAR-2.7.11b

# Get human genome files and start building indices
# Because read length could vary, setting the sjdbOverhang setting to 100bp
# Note that if running this locally (won't be necessary if running on a VM), you may need to change the 
# docker disk location so that you have enough disk space to build the genome
# To do this locally (on WSL2, make sure that you have this set up): 
#       Backup/move files on NTFS USB SSD
#       In PowerShell (PS), use `wmic diskdrive list brief` to ID the DeviceID of the USB SSD
#       (PS) `wsl --mount \\.\PHYSICALDRIVE` to attach the device
#       In WSL (WSL), use `sudo fdisk -l` to ID where the USB SSD is attached to (should be /dev/sdX)
#       (WSL) Confirm location by `lsblk`, then format the drive as ext4 using `sudo mkfs.ext4 /dev/sdX`
#       (WSL) Create and mount the drive using `sudo mkdir -p /mnt/ext_drive && sudo mount /dev/sdX /mnt/ext_drive`
#       (WSL) Verify it's mounted using `df -h /mnt/ext_drive`
#       (WSL) Now migrate docker over to new location using `mkdir -p /mnt/ext_drive/docker && \
#                                                            sudo cp -a /var/lib/docker/ /mnt/ext_drive/docker/`
#       (WSL) If a daemon file doesn't exist (`less /etc/docker/daemon.json`), create one using
#             `echo { "data-root": "/mnt/ext_drive/docker" } >> daemon.json`, then migrate that file into 
#             /etc/docker by `sudo mv daemon.json /etc/docker/`
#       (WSL) Don't forget to cleanup data (if desired) from /var/lib/docker/
# Now when you restart Docker, it should be running out of the USB SSD and have your images/containers ready

RUN cd /home/apps/STAR/index && \
    wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz && \
    wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz && \
    gunzip Homo_sapiens.GRCh38.113.gtf.gz Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# Running genome generation step with less RAM, change this depending on the system configuration
RUN STAR --runMode genomeGenerate --runThreadN 8 --genomeDir /home/apps/STAR/index --genomeFastaFiles /home/apps/STAR/index/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --sjdbGTFfile /home/apps/STAR/index/Homo_sapiens.GRCh38.113.gtf --sjdbOverhang 100 --limitGenomeGenerateRAM 14000000000 --genomeChrBinNbits 16 --genomeSAsparseD 2 

# Cleanup genome generation files
RUN rm /home/apps/STAR/index/Homo_sapiens.GRCh38.113.gtf /home/apps/STAR/index/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Set environment variables for STAR
RUN echo "STAR_GENOME=/home/apps/STAR/index" >> /etc/environment && \
    source /etc/environment