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

// Check that NF_WORKDIR is defined in environment
def nf_workdir = System.getenv('NF_WORKDIR')
if (!nf_workdir) {
    log.error "Environment variable NF_WORKDIR is not set."
    log.info  "Please set it, e.g.: export NF_WORKDIR=/path/to/atlas"
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
def results_dir = file("${nf_core_bulk_quantification}/${params.EXP_ID}")
results_dir.mkdirs()



// -------------------- PROCESSES -------------------- //

// include a process that checks goofys mount, mounts if non-existent

process GET_SAMPLES {
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


// Get species information
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

// Run the nf-core/rnasesq workflow
process RUN_RNASEQ {

    publishDir params.outdir, mode: 'copy'

    input:
    path samplesheet
    val  EXP_ID
    path params_json, stageAs: "${EXP_ID}_params.json"

    output:
    path "${EXP_ID}.rnaseq.done"

    script:
    """
    set -euo pipefail

    echo "Running RNA-seq subworkflow for ${EXP_ID}"

    # Extract FASTA path from JSON
    if command -v jq &> /dev/null; then
        FASTA_PATH=\$(jq -r '.fasta // .genome // empty' "${params_json}")
    else
        FASTA_PATH=\$(grep -oP '"fasta"\\s*:\\s*"\\K[^"]+' "${params_json}" || \
                     grep -oP '"genome"\\s*:\\s*"\\K[^"]+' "${params_json}")
    fi
    
    GENOME_FASTA_INDEX="\${FASTA_PATH}.fai"
    
    # Check if CSI is needed
    if [ -f "\$GENOME_FASTA_INDEX" ]; then
        if awk '\$2 > 512000000 {exit 1}' "\$GENOME_FASTA_INDEX"; then
            BAM_INDEX=""
            echo "BAI index (chromosomes <512 Mbp)"
        else
            BAM_INDEX="--bam_csi_index"
            echo "CSI index (chromosomes >512 Mbp detected)"
        fi
    else
        BAM_INDEX=""
        echo "FASTA index not found: \$GENOME_FASTA_INDEX (defaulting to BAI)"
    fi

    nextflow run ${workflow.projectDir}/subworkflows/rnaseq/main.nf \\
        -params-file "${params_json}" \\
        -c "${workflow.projectDir}/conf/rnaseq.config" \\
        -profile singularity \\
        \$BAM_INDEX \\
        --without-wave \\
        --skip_fastqc  \\
        --skip_rseqc  \\
        --skip_qualimap \\
        --skip_dupradar \\
        --skip_preseq \\
        --skip_biotype_qc \\
        --skip_kraken2 \\
        --skip_stringtie \\
        --skip_deseq2_qc \\
        --skip_markduplicates \\
        --skip_bigwig \\
        -with-trace "${params.outdir}/${EXP_ID}_trace.tsv" \\
    && nextflow clean -f

    # Create done file only if workflow succeeded
    touch "${EXP_ID}.rnaseq.done"
    
    """
}

process MULTIQC_SANITISATION {

    publishDir params.outdir, mode: 'copy'

    output:
    path "multiqc_sanitisation.done"

    script:
    """
    set -euo pipefail

    REPORT_DIR="${params.outdir}/multiqc/star_salmon"
    INPUT_HTML="\$REPORT_DIR/multiqc_report.html"
    BACKUP_HTML="\$REPORT_DIR/multiqc_report_original.html"
    OUTPUT_HTML="\$REPORT_DIR/multiqc_report.html"

    cp "\$INPUT_HTML" "\$BACKUP_HTML"

    sed \\
      \${REF_ESC:+-e "s#\${referencePath}#REFERENCES_PATH#g"} \\
      \${WORK_ESC:+-e "s#\${nf_workdir}#WORKDIR#g"} \\
      \${OUT_ESC:+-e "s#\${params.outdir}#OUT_DIR#g"} \\
      \${PROJ_ESC:+-e "s#\${workflow.projectDir}#GIT-REPO#g"} \\
      "\$BACKUP_HTML" > "\$OUTPUT_HTML"

    touch multiqc_sanitisation.done
    """
}

// For cleanup
process HANDLE_STATUS {

    publishDir params.outdir, mode: 'copy'
    
    input:
    val pipeline_status
    
    script:
    """
    echo "Cleaning up for: ${params.EXP_ID}"

    if [ "${pipeline_status}" == "SUCCESS" ]
        # Success flag for downstream logic / idempotency
        echo "Creating ${params.EXP_ID}.rnaseq.done"
        touch "${params.EXP_ID}.rnaseq.done"
    else
        # Mark failure
        echo "Creating ${params.EXP_ID}.rnaseq.fail"
        touch "${params.EXP_ID}.rnaseq.fail"
        # Grab error from log, write to file
        errOut=\$( echo "Unknown error" ) ### Command here to grab error ################
        echo \$errOut >> ${nf_core_bulk_quantification}/excluded.txt
    fi

    """
}


workflow {
    samplesheet = GET_SAMPLES(params.EXP_ID)
    params_json_ch = GET_SPECIES(params.EXP_ID)
    rnaseq_done = RUN_RNASEQ(samplesheet, params.EXP_ID, params_json_ch)
    MULTIQC_SANITISATION(rnaseq_done)
}

workflow sanitise_only {
    MULTIQC_SANITISATION()
}

workflow.onComplete {
    log.info "Pipeline completed at $workflow.complete."

    if (workflow.success) {
        log.info "Pipeline completed successfully!"
    } else {
        log.info "Pipeline failed with exit status: ${workflow.exitStatus}"
        def err_report = workflow.errorReport?.toString()
        def err_msg = workflow.errorMessage
        def failureFile = file("${params.outdir}/${params.EXP_ID}.rnaseq.failed")
        failureFile.text = "${err_msg} \n\n ${err_report}"

        // Write to excluded.txt
        def excluded = file("${nf_core_bulk_quantification}/excluded.txt")   
        excluded << "${params.EXP_ID}\t${failureFile}\n"
    }
}

