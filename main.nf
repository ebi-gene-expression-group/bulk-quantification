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


// Set experiment specific parameters
process SET_PARAMS {

    input:
        val EXP_ID

    output:
        path "${EXP_ID}_params.json"

    script:
    """
    bash ${projectDir}/bin/generate_params.sh ${EXP_ID}
    """
}

// Get species information
process GET_TAX_ID {

    input:
        val EXP_ID

    output:
        stdout emit: TAX_ID

    script:
    """
    bash ${projectDir}/bin/generate_params.sh ${EXP_ID}
    """
}

// Run the nf-core/rnasesq workflow
process RUN_RNASEQ {

    publishDir params.outdir, mode: 'copy'

    input:
    path samplesheet
    val  EXP_ID
    val TAX_ID
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

    if command -v jq &> /dev/null; then
        RIBO_INDEX=\$(jq -r '.ribo_database_index // empty' "${params_json}")
        RIBO_MANIFEST=\$(jq -r '.ribo_database_manifest // empty' "${params_json}")
        CONTAM_INDEX=\$(jq -r '.contamination_index // empty' "${params_json}")
    else
        RIBO_INDEX=\$(grep -oP '"ribo_database_index"\\s*:\\s*"\\K[^"]+' "${params_json}" || true)
        RIBO_MANIFEST=\$(grep -oP '"ribo_database_manifest"\\s*:\\s*"\\K[^"]+' "${params_json}" || true)
        CONTAM_INDEX=\$(grep -oP '"contamination_index"\\s*:\\s*"\\K[^"]+' "${params_json}" || true)
    fi

    if [[ -z "\${RIBO_INDEX}" ]]; then
        echo "Missing required ribo_database_index in ${params_json}"
        exit 1
    fi

    if [[ -z "\${RIBO_MANIFEST}" ]]; then
        echo "Missing required ribo_database_manifest in ${params_json}"
        exit 1
    fi

    if [[ -z "\${CONTAM_INDEX}" ]]; then
        echo "Missing required contamination_index in ${params_json}"
        exit 1
    fi
    
    # Check if ribo database index directory exists AND is not empty
    if [ -d "\${RIBO_INDEX}" ] && [ "\$(ls -A "\${RIBO_INDEX}")" ]; then
        echo "Using existing SortMeRNA index from: \${RIBO_INDEX}"
    else
        echo "SortMeRNA index not found. Create a new index..."
        exit 1
    fi

    if [ -f "\${RIBO_MANIFEST}" ]; then
        echo "Using existing SortMeRNA manifest from: \${RIBO_MANIFEST}"
        cat \${RIBO_MANIFEST}
        missing=0
        while IFS= read -r f; do
            [[ -z "\$f" ]] && continue
            if [[ ! -e "\$f" ]]; then
                echo "Missing: \$f"
                missing=1
            fi
        done < "\${RIBO_MANIFEST}"
    
        if [[ \$missing -eq 0 ]]; then
            echo "All files exist."
        else
            echo "Some files missing and sortmerna likely to fail, exiting..."
            exit 1
        fi
    else
        echo "SortMeRNA manifest not found. Create a new manifest..."
        exit 1
    fi

    # Check if contamination index directory exists AND is not empty
    if [ -d "\${CONTAM_INDEX}" ] && [ "\$(ls -A "\${CONTAM_INDEX}")" ]; then
        echo "Using existing contamination index from: \${CONTAM_INDEX}"
    else
        echo "Contamination index directory not found or empty: \${CONTAM_INDEX}"
        exit 1
    fi

    # Get STAR profile name from Python script
    STAR_PROFILE=\$(python "${workflow.projectDir}/bin/tax_id_to_profile.py" "${TAX_ID}" | tr -d '[:space:]')

    if [[ -z "\${STAR_PROFILE}" ]]; then
        echo "Python script did not return a STAR profile, using default"
        STAR_PROFILE="default"
    fi

    echo "Using STAR profile: \${STAR_PROFILE}"

    nextflow run ${workflow.projectDir}/subworkflows/rnaseq/main.nf \\
        -params-file "${params_json}" \\
        -c "${workflow.projectDir}/conf/rnaseq.config" \\
        -c "${workflow.projectDir}/conf/star_\${STAR_PROFILE}.config" \\
        -profile singularity \\
        \$BAM_INDEX \\
        --contaminant_screening kraken2_bracken \\
        --kraken_db "\${CONTAM_INDEX}" \\
        --without-wave \\
        --save_unaligned \\
        --skip_bbsplit \\
        --skip_fastqc \\
        --skip_rseqc \\
        --skip_qualimap \\
        --skip_dupradar \\
        --skip_preseq \\
        --skip_biotype_qc \\
        --skip_stringtie \\
        --skip_deseq2_qc \\
        --skip_markduplicates \\
        --skip_bigwig \\
        --remove_ribo_rna \\
        --ribo_removal_tool sortmerna \\
        --ribo_database_manifest "\${RIBO_MANIFEST}" \\
        --sortmerna_index "\${RIBO_INDEX}" \\
        -with-trace "${params.outdir}/${EXP_ID}_trace.tsv" \\
        -with-tower \\
        -name "nf_core_rnaseq_${EXP_ID}"

# \\
#    && nextflow clean -f
    
    # Create done file only if workflow succeeded
    touch "${EXP_ID}.rnaseq.done"
    
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
    params_json_ch = SET_PARAMS(params.EXP_ID)
    GET_TAX_ID(params.EXP_ID)
    TAX_ID = GET_TAX_ID.out.TAX_ID.map { it.trim() }
    RUN_RNASEQ(samplesheet, params.EXP_ID, TAX_ID, params_json_ch)
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
