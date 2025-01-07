#!/usr/bin/env nextflow

/* 
Nextflow script to run processing on bulk RNA-seq data
    using Docker container previously created (see build_rnaseq_environ.nf).
    The docker container will be named rnaseq_environ.

Adapted from the Cebola Lab's RNA-seq pipeline posted here:
    https://github.com/CebolaLab/RNA-seq

Some test data to use can be found at:https://www.ebi.ac.uk/ena/browser/view/PRJEB34752
    and is part of the study found at: https://www.nature.com/articles/s41588-022-01100-4
    "Single-cell and bulk transcriptome sequencing identifies two epithelial tumor cell 
    states and refines the consensus molecular classification of colorectal cancer"
    by Juanito I et al (2022). The paper will include patient details for the samples.
The data includes both border (normal) and core (tumor) samples from colorectal cancer patients.
It is single-end, 51bp reads.
*/

params.input = "/path/to/fastq_files/*.fastq.gz"    // Raw data direectory
//params.genome = "/path/to/genome/index"              // Genome index directory inside rnaseq_environ
params.output = "./output"                         // Output directory

log.info """\
    -----------------------------------------
    Bulk RNA-seq Data : Processing (Nextflow)
    -----------------------------------------
    Processing bulk RNA-seq data using Docker container
    (rnaseq_environ) previously created. Files should be
    in fastq.gz format.

    input:          $params.input
    genome index:   $params.genome
    output:         $params.output
    """
    .stripIndent(true)

// Define the workflow
workflow {
    processFastqFiles(params.input, params.genome, params.output)
}

// Define the process
process processFastqFiles {
    container 'my_docker_image:latest'          // Use the prebuilt Docker container

    input:
    path fastq_files from Channel.fromPath(params.input) // Input FASTQ files

    output:
    path 'output/*' into processed_results       // Output results directory

    script:
    """
    mkdir -p output
    my_tool --input ${fastq_files} \               # Command inside the container
            --index /home/apps/STAR/index \        # Genome index location in rnaseq_environ
            --output output/
    """
}

// Save processed results to a local directory
processed_results.view { file -> "Processed file: ${file}" }