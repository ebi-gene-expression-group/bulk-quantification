#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// -------------------- PARAM CHECKS -------------------- //

// Require EXP_ID
if (!params.EXP_ID) {
    log.error "Missing required parameter: --EXP_ID"
    log.info  "Usage: nextflow run main.nf --EXP_ID <experiment_id>"
    System.exit(1)
}

// Check that ATLAS_PROD is defined in environment
def atlasProd = System.getenv('ATLAS_PROD')
if (!atlasProd) {
    log.error "Environment variable ATLAS_PROD is not set."
    log.info  "Please set it, e.g.: export ATLAS_PROD=/path/to/atlas"
    System.exit(1)
}

def nf_core_bulk_quantification = System.getenv('NF_CORE_BULK_QUANTIFICATION')

if (!nf_core_bulk_quantification) {
    log.error "Environment variable NF_CORE_BULK_QUANTIFICATION is not set."
    log.info  "Please set it, e.g.: export NF_CORE_BULK_QUANTIFICATION=/path/to/atlas"
    System.exit(1)
}


// Check that REFERENCES_PATH is defined in environment
def referencePath = System.getenv('BULK_REFERENCES_DIR')
if (!referencePath) {
    log.error "Environment variable BULK_REFERENCES_DIR is not set."
    log.info  "Please set it, e.g.: export REFERENCES_PATH=/path/to/atlas"
    System.exit(1)
}


def era_public_mount_path = System.getenv('ERA_PUBLIC_MOUNT_PATH')
if (!era_public_mount_path) {
    log.error "Environment variable ERA_PUBLIC_MOUNT_PATH  is not set."
    log.info  "Please set it, e.g.: export ERA_PUBLIC_MOUNT_PATH=/path/to/atlas"
    System.exit(1)
}

def fastq_rawdata_dir = System.getenv('FASTQ_RAWDATA_DIR')
if (!fastq_rawdata_dir) {
    log.error "Environment variable FASTQ_RAWDATA_DIR  is not set."
    log.info  "Please set it, e.g.: export FASTQ_RAWDATA_DIR=/path/to/atlas"
    System.exit(1)
}

// Consider a process logic that does require defining the variables below
def endpoint_url = System.getenv('FIRE_ENDPOINT')
if (!fastq_rawdata_dir) {
    log.error "Environment variable FIRE_ENDPOINT is not set."
    log.info  "Please set it, e.g.: export FIRE_ENDPOINT=https://<hostname>/"
    System.exit(1)
}

def era_public_s3_path = System.getenv('ERA_PUBLIC_S3_PATH')
if (!fastq_rawdata_dir) {
    log.error "Environment variable ERA_PUBLIC_S3_PATH is not set."
    log.info  "Please set it, e.g.: export ERA_PUBLIC_S3_PATH=s3://path/dir"
    System.exit(1)
}

// Define output directory based on EXP_ID
params.outdir = "${nf_core_bulk_quantification}/${params.EXP_ID}"

workflow {
    samplesheet = get_samples(params.EXP_ID)
    params_json_ch = GET_SPECIES(params.EXP_ID)
    run_rnaseq(samplesheet, params.EXP_ID, params_json_ch)
}


// -------------------- PROCESSES -------------------- //

// include a process that checks goofys mount, mounts if non-existent

process get_samples {
    //container "$params.aws_container"

    input:
    val EXP_ID

    output:
    path "${EXP_ID}_samplesheet.csv"

    script:
    """
    export EXP_ID=${EXP_ID}
    export CONFIG="atlas"
    export SAMPLESHEET="\${EXP_ID}_samplesheet.csv"

    echo 'Creating samplesheet for \${EXP_ID}'

    bash '${projectDir}/bin/get_samples.sh' \\
        -a "\${EXP_ID}" \\
        -x "\${CONFIG}" \\
        -s "\${SAMPLESHEET}" \\
        -e "${endpoint_url}" \\
        -m "${era_public_s3_path}" \\
        -c "${fastq_rawdata_dir}"
    """
}


process GET_SPECIES {
  input:
    val EXP_ID

  output:
    path "${EXP_ID}_params.json"

  script:
  """
  ${projectDir}/bin/generate_params.sh ${EXP_ID}
  """
}


process run_rnaseq {
    publishDir params.outdir, mode: 'copy'

    input:
    path samplesheet
    val  EXP_ID
    path "${EXP_ID}_params.json"

    output:
    path "${EXP_ID}.rnaseq.done"

    script:
    """
    set -euo pipefail

    echo "Running RNA-seq subworkflow for \${EXP_ID}"

    nextflow run ${projectDir}/subworkflows/rnaseq/main.nf \\
        -params-file "\${EXP_ID}_params.json" \\
        -c "${projectDir}/conf/rnaseq.config" \\
        -profile singularity \\
        --without-wave \\
        -with-trace "${params.outdir}/${EXP_ID}_trace.tsv"

    exit_code=\$?

    if [ \$exit_code -eq 0 ]; then
        # Success flag for downstream logic / idempotency
        touch "\${EXP_ID}.rnaseq.done"
    else
        # Mark failure
        touch "\${EXP_ID}.rnaseq.fail"
        # Grab error from log, write to file
        errOut=\$( echo "Unknown error" ) ### Command here to grab error ################
        echo \$errOut >> ${nf_core_bulk_quantification}/excluded.txt
    fi
    """
}
