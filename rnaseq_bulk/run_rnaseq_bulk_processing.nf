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

params.reads = "/path/to/fastq_files/*.fastq.gz"    // Raw data directory
params.docker = "rnaseq_environ"                    // Docker container name  
params.genome = "/home/apps/STAR/index"             // Genome index directory inside docker container
params.outdir = "/path/to/output"                   // Output directory

log.info """\
    -----------------------------------------
    Bulk RNA-seq Data : Processing (Nextflow)
    -----------------------------------------
    Processing bulk RNA-seq data using Docker container
    (rnaseq_environ) previously created. Files should be
    in fastq.gz format.

    input:              $params.reads
    docker container:   $params.docker
    genome index:       $params.genome
    output:             $params.outdir
    """
    .stripIndent(true)

Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .set{ input_ch }

input_ch.view()

// Workflow <fill>
/*
workflow {
    processFastqFiles(params.input, params.genome, params.docker, params.output)
}

// Create a channel of read files
process processFastqFiles {
    container params.docker          // Use the prebuilt Docker container

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
*/

#!/usr/bin/env nextflow

params.reads = '/home/byuen/projects/bioinformatics_projects/rnaseq_bulk/test_data/*_r*.fastq.gz'   // Raw data directory
params.docker = 'rnaseq_environ'                                                                    // Docker container name  
params.genome = '/home/apps/STAR/index'                                                             // Genome index directory inside docker container
params.outdir = 'output'                                                                            // Output directory

log.info """\
    -----------------------------------------
    Bulk RNA-seq Data : Processing (Nextflow)
    -----------------------------------------
    Processing bulk RNA-seq data using Docker container
    (rnaseq_environ) previously created. Files should be
    in fastq.gz format.

    input:              $params.reads
    docker container:   $params.docker
    genome index:       $params.genome
    output:             $params.outdir
    """
    .stripIndent(true)

println "Starting"

nextflow.enable.dsl=2

params.input = '/home/byuen/projects/bioinformatics_projects/rnaseq_bulk/test_data/*_r*.fastq.gz'

process readFastqFiles {
    input:
    tuple val(sample_id), path(reads) from collectFiles(params.input)

    output:
    path("${sample_id}_output.txt") into results

    script:
    if (reads.size() == 1) {
        """
        echo "Single-end sample: $sample_id with file ${reads[0].getName()}" > ${sample_id}_output.txt
        """
    } else if (reads.size() == 2) {
        """
        echo "Paired-end sample: $sample_id with files ${reads[0].getName()} and ${reads[1].getName()}" > ${sample_id}_output.txt
        """
    }
}

workflow {
    results = collectFiles(params.input)
        | readFastqFiles

    results.view()
}

// Helper function to collect files into sets
def collectFiles(pattern) {
    Channel
        .fromFilePairs(pattern, flat: false)
        .map { id, files -> [id, files] }
}


/*workflow {
    files_ch = test(reads_ch)
    files_ch.view()
}

process test {
    input:
    val reads

    output:
    stdout

    script:
    """
    echo "testing: $reads"
    """
}
*/