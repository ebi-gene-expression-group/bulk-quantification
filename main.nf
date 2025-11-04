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

// Check that REFERENCES_PATH is defined in environment
def referencePath = System.getenv('REFERENCES_PATH')
if (!referencePath) {
    log.error "Environment variable REFERENCES_PATH is not set."
    log.info  "Please set it, e.g.: export REFERENCES_PATH=/path/to/atlas"
    System.exit(1)
}

// Define output directory based on EXP_ID and ATLAS_PROD
params.outdir = "${atlasProd}/analysis/baseline/rnaseq/experiments/${params.EXP_ID}"

workflow {
    samplesheet = create_samplesheet(params.EXP_ID)
    // run_rnaseq(samplesheet, params.EXP_ID)
}


// -------------------- PROCESSES -------------------- //

process create_samplesheet {
    publishDir "${params.outdir}/samplesheet", mode: 'copy'

    input:
    val EXP_ID

    output:
    path "samplesheet.csv"

    script:
    """
    echo "Creating samplesheet for ${EXP_ID}"
    grep "<assay>" ${params.outdir}/${EXP_ID}-configuration.xml | sed 's/\s*<\/*assay>//g' > ${EXP_ID}_ids.csv
    /bin/create_samplesheet.sh ${EXP_ID}_ids.csv > ${EXP_ID}_samplesheet.csv
    """
}


process run_rnaseq {
    publishDir "${params.outdir}/rnaseq", mode: 'copy'

    input:
    path samplesheet
    val  EXP_ID

    output:
    path "${EXP_ID}.rnaseq.flag"

    script:
    """
    set -euo pipefail

    # Prepare params file dynamically
    export EXP_ID="${EXP_ID}"

    echo "Running RNA-seq subworkflow for ${EXP_ID}"

    envsubst < params.template.json > ${EXP_ID}_params.json

    nextflow run subworkflows/rnaseq/main.nf \\
        -params-file ${EXP_ID}_params.json \\
        -C conf/rnaseq.config \\
        -profile singularity \\
        --without-wave

    # Mark success
    touch ${EXP_ID}.rnaseq.done
    """
}
