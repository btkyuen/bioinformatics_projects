#!/usr/bin/env nextflow

process sayHello {
  input: 
    val x
  output:
    stdout
  script:
    """
    docker run -di --name image_rnaseq_bulk ubuntu:20.04 bin/bash && \
    docker exec image_rnaseq_bulk bash -c "apt-get update && \
        DEBIAN_FRONTEND=noninteractive && \
        apt-get install -y -q curl \
            wget \
            gfortran \
            build-essential \
            tzdata \
            software-properties-common && \
        wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
        add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
        apt update && \
        apt-get install -y --no-install-recommends \
            r-base-core=4.3.3-1.2004.0 \
            r-base-html=4.3.3-1.2004.0 \
            r-doc-html=4.3.3-1.2004.0 \
            r-base-dev=4.3.3-1.2004.0"
    """
}

workflow {
  Channel.of('Bonjour', 'Ciao', 'Hello', 'Hola') | sayHello | view
}