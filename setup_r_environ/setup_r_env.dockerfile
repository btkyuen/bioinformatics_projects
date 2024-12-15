# Use Ubtunu 20.04 and R 4.3.3
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Set up Ubuntu requirements for R and R packages
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
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libudunits2-dev \
        libgdal-dev \
        libgeos-dev \
        libproj-dev \
        libmagick++-dev \
        libavfilter-dev \
        libpoppler-cpp-dev \
        libtesseract-dev \
        libleptonica-dev \
        cargo \
        tesseract-ocr-eng

# Install R including key and repo
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
    apt update && \
    apt-get install -y -q --no-install-recommends \
        r-base-core=4.3.3-1.2004.0 \
        r-base-html=4.3.3-1.2004.0 \
        r-doc-html=4.3.3-1.2004.0 \
        r-base-dev=4.3.3-1.2004.0 && \
    apt-get clean

# Install R packages and dependencies
RUN Rscript -e "install.packages(c('dplyr', 'textshaping', 'ragg'), dependencies = TRUE)"
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-5.tar.gz', dependencies = TRUE)"
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-60.tar.gz', dependencies = TRUE)"
RUN Rscript -e "install.packages(c('lattice', 'reticulate', 'sf', 'ggplot2', 'tidyplot', 'stringr', 'future', 'future.apply', \
                                    'progressr', 'spam'), dependencies = TRUE)"
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_5.0.1.tar.gz', dependencies = TRUE)"
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/fitdistrplus/fitdistrplus_1.1-11.tar.gz', dependencies = TRUE)"
RUN Rscript -e "install.packages(c('leiden', 'av', 'pdftools', 'gifski', 'tesseract', 'Cairo', 'cowplot', 'fastDummies', 'ggrepel', 'ggridges', \
                    'ica', 'irlba', 'lmtest', 'matrixStats', 'plotly', 'RANN', 'RcppAnnoy', 'RcppHNSW', 'ROCR', 'RSpectra', \
                    'Rtsne', 'uwot', 'sctransform', 'scattermore'), dependencies = TRUE)"
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_5.0.3.tar.gz', dependencies = TRUE)"
