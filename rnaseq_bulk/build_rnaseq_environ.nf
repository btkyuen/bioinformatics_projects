#!/usr/bin/env nextflow

/* 
Nextflow script to build a Docker image for RNA-seq analysis
and start the container for next steps in the data processing.
*/

log.info """\
    ------------------------------------------
    Bulk RNA-seq Data : Environment (Nextflow)
    ------------------------------------------
    Building Docker container for RNA-seq analysis.
    container:    rnaseq_environ
    """
    .stripIndent(true)

workflow{
  buildRNAseqEnviron()
  run_buildRNAseqEnviron()
}

process buildRNAseqEnviron{
  input:
  path 'Dockerfile'

  output:
  file 'environ_build.log'

  script:
  """
  docker build -f=setup_rnaseq_environ.dockerfile -t=rnaseq_environ . > environ_build.log 2>&1
  echo "Docker image built successfully." >> environ_build.log
  """
}

process run_buildRNAseqEnviron{
  input:
  path 'environ_build.log'

  output:
  file 'environ_run.log'

  script:
  """
  docker run rnaseq_environ:latest > environ_run.log 2>&1
  echo "Container run completed." >> environ_run.log
  """
}