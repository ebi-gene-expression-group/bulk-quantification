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


// Check that BULK_REFERENCES_DIR is defined in environment
def referencePath = System.getenv('BULK_REFERENCES_DIR')
if (!referencePath) {
    log.error "Environment variable BULK_REFERENCES_DIR is not set."
    log.info  "Please set it, e.g.: export BULK_REFERENCES_DIR=/path/to/references"
    System.exit(1)
}

// Check that NF_WORKDIR is defined in environment
def nf_workdir = System.getenv('NF_WORKDIR')
if (!nf_workdir) {
    log.error "Environment variable NF_WORKDIR is not set."
    log.info  "Please set it, e.g.: export NF_WORKDIR=/path/to/workdir"
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
if (!endpoint_url) {
    log.error "Environment variable FIRE_ENDPOINT is not set."
    log.info  "Please set it, e.g.: export FIRE_ENDPOINT=https://<hostname>/"
    System.exit(1)
}

def era_public_s3_path = System.getenv('ERA_PUBLIC_S3_PATH')
if (!era_public_s3_path) {
    log.error "Environment variable ERA_PUBLIC_S3_PATH is not set."
    log.info  "Please set it, e.g.: export ERA_PUBLIC_S3_PATH=s3://path/dir"
    System.exit(1)
}

// Define output directory based on EXP_ID
params.outdir = "${nf_core_bulk_quantification}/${params.EXP_ID}"
def results_dir = file(params.outdir)
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
        path "${EXP_ID}_params.json", emit: params_json
        path "tax_id.txt",            emit: tax_id

    script:
    """
    bash ${projectDir}/bin/generate_params.sh ${EXP_ID} > tax_id.txt
    """
}

// Get species information
process GET_STAR_PROFILE {

    container "quay.io/ebigxa/ete3:v1.0"

    containerOptions = "--bind ${referencePath}:${referencePath}"

    input:
        path tax_id_file

    output:
        path "*.config"

    script:
    """
    set -euo pipefail

    ########## TESTS ##########
    echo "BELOW ARE MOUNTED PATH TESTS"
    ls -l "${referencePath}/taxonomy/taxa.sqlite"
    ls -l "${workflow.projectDir}" | head
    echo "MOUNTED PATH TESTS DONE"
    ########## TESTS ##########

    TAX_ID=\$(cat "${tax_id_file}" | tr -d '[:space:]')

    # Get STAR profile name from Python script; non-zero exit aborts the process
    STAR_PROFILE=\$(python "${workflow.projectDir}/bin/tax_id_to_profile.py" "\${TAX_ID}" | tr -d '[:space:]')

    STAR_CONFIG="${workflow.projectDir}/conf/star_\${STAR_PROFILE}.config"

    echo "Using STAR profile: \$(basename "\${STAR_CONFIG}")"

    if [[ ! -f "\${STAR_CONFIG}" ]]; then
        echo "ERROR: STAR profile config not found: \${STAR_CONFIG}" >&2
        exit 1
    fi

    cp "\${STAR_CONFIG}" .
    """
}

// Run the nf-core/rnasesq workflow
process RUN_RNASEQ {

    publishDir params.outdir, mode: 'copy'

    input:
    path samplesheet
    val  EXP_ID
    path star_config
    path params_json, stageAs: "${EXP_ID}_params.json"

    output:
    path "${EXP_ID}.rnaseq.done"
    path "star_profile.log"
    path "multiqc/star_salmon/multiqc_report.html"

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

    STAR_PROFILE_USED="\$(basename "${star_config}")"
    mkdir -p "${params.outdir}"
    printf 'STAR_PROFILE_USED\t%s\n' "\${STAR_PROFILE_USED}" > "${params.outdir}/star_profile.log"
    cp "${params.outdir}/star_profile.log" star_profile.log

    if [ ! -f "${workflow.projectDir}/subworkflows/rnaseq/main.nf" ]; then
        echo "ERROR: rnaseq subworkflow not found at ${workflow.projectDir}/subworkflows/rnaseq/main.nf. Did you run 'git submodule update --init --recursive'?" >&2
        exit 1
    fi

    nextflow run ${workflow.projectDir}/subworkflows/rnaseq/main.nf \\
        -params-file "${params_json}" \\
        -c "${workflow.projectDir}/conf/rnaseq.config" \\
        -c "${star_config}" \\
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
        --multiqc_config ${workflow.projectDir}/conf/multiqc_star_profile.yaml \\
        -with-trace "${params.outdir}/${EXP_ID}_trace.tsv" \\
        -with-tower \\
        -name "nf_core_rnaseq_${EXP_ID}"

    run_status=\$?
    if [[ \$run_status -eq 0 ]]; then
        nextflow clean -f
    fi
    
    # Stage multiqc HTML into task work dir for downstream processes
    MULTIQC_HTML="${params.outdir}/multiqc/star_salmon/multiqc_report.html"
    if [ ! -f "\$MULTIQC_HTML" ]; then
        echo "ERROR: MultiQC report not found at \$MULTIQC_HTML" >&2
        exit 1
    fi
    mkdir -p multiqc/star_salmon
    cp "\$MULTIQC_HTML" multiqc/star_salmon/multiqc_report.html

    # Create done file only if workflow succeeded
    touch "${EXP_ID}.rnaseq.done"
    
    """
}

process MULTIQC_SANITISATION {

    publishDir "${params.outdir}/multiqc/star_salmon", mode: 'copy', overwrite: true, pattern: "multiqc_report*.html"
    publishDir params.outdir, mode: 'copy', pattern: "multiqc_sanitisation.done"

    input:
    path multiqc_html

    output:
    path "multiqc_report.html"
    path "multiqc_report_original.html"
    path "multiqc_sanitisation.done"

    script:
    // Build sed expressions in Groovy
    def esc = { it.toString().replaceAll(/([\\#&])/,'\\\\$1') }

    def sed_cmds = []
    
    if (referencePath)
        sed_cmds << "-e 's#${esc(referencePath)}#BULK_REFERENCES_DIR#g'"
    
    if (nf_workdir)
        sed_cmds << "-e 's#${esc(nf_workdir)}#WORKDIR#g'"
    
    if (params.outdir)
        sed_cmds << "-e 's#${esc(params.outdir)}#OUT_DIR#g'"
    
    if (workflow.projectDir)
        sed_cmds << "-e 's#${esc(workflow.projectDir)}#GIT-REPO#g'"
    
    def sed_string = sed_cmds.join(' ')

    """
    set -euo pipefail

    cp "${multiqc_html}" multiqc_report_original.html
    sed ${sed_string} multiqc_report_original.html > multiqc_report.html

    touch multiqc_sanitisation.done
    """
}

workflow {
    samplesheet = GET_SAMPLES(params.EXP_ID)
    SET_PARAMS(params.EXP_ID)
    star_config_ch = GET_STAR_PROFILE(SET_PARAMS.out.tax_id)
    params_json_ch = SET_PARAMS.out.params_json
    (rnaseq_done, multiqc_html) = RUN_RNASEQ(samplesheet, params.EXP_ID, star_config_ch, params_json_ch)
    MULTIQC_SANITISATION(multiqc_html)
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
