#!/bin/bash

mkdir -p /home/references/STAR/index/ && \
    mkdir -p /home/references/salmon/index/

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

cd /home/references && \
    wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz && \
    wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz && \
    gunzip Homo_sapiens.GRCh38.113.gtf.gz Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

STAR --runMode genomeGenerate \
        --runThreadN 8 \
        --genomeDir /home/references/STAR/index \
        --genomeFastaFiles /home/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        --sjdbGTFfile /home/references/Homo_sapiens.GRCh38.113.gtf \
        --sjdbOverhang 100 \
        --limitGenomeGenerateRAM 14000000000 \
        --genomeChrBinNbits 16 \
        --genomeSAsparseD 2 

# Generate Salmon index
gffread -w /home/references/salmon/index/Homo_sapiens.GRCh38.dna.primary_assembly.transcripts.fa \
        -g /home/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        /home/references/Homo_sapiens.GRCh38.113.gtf && \
    mv /home/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai /home/references/salmon/index/

# Set environment variables for Salmon and STAR
echo 'SALMON_GENOME=/home/references/salmon/index/Homo_sapiens.GRCh38.dna.primary_assembly.transcripts.fa' >> /etc/environment && \
    echo 'STAR_GENOME=/home/references/STAR/index' >> /etc/environment && \
    source /etc/environment

# Cleanup genome generation files
rm /home/references/Homo_sapiens.GRCh38.113.gtf /home/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa

