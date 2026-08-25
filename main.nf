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

params.batch_size = (params.batch_size ?: 0) as Integer
if (params.batch_size < 0) {
    log.error "Invalid --batch_size '${params.batch_size}'. Use 0 to disable batching, or a positive integer."
    System.exit(1)
}

// Resume process RUN_RNASEQ?
def resumeOpt = System.getenv('RESUME_RNASEQ') ?: ""


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

process CREATE_BATCH_SAMPLESHEETS {

    input:
    path samplesheet
    val  EXP_ID
    val  batch_size

    output:
    path "batch_manifest.tsv", emit: manifest
    path "batch_samplesheets/*.csv", emit: sheets

    script:
    """
    set -euo pipefail

    mkdir -p batch_samplesheets

    header=\$(head -n1 "${samplesheet}")
    total_lines=\$(wc -l < "${samplesheet}")
    total_samples=\$(( total_lines - 1 ))

    if [[ \$total_samples -le 0 ]]; then
        echo "ERROR: Empty samplesheet generated for ${EXP_ID}: ${samplesheet}" >&2
        exit 1
    fi

    if [[ ${batch_size} -le 0 || \$total_samples -le ${batch_size} ]]; then
        batch_csv="batch_samplesheets/${EXP_ID}_batch_000001_samplesheet.csv"
        cp "${samplesheet}" "\$batch_csv"
    else
        tail -n +2 "${samplesheet}" | split -l ${batch_size} -d -a 6 --numeric-suffixes=1 - "batch_samplesheets/${EXP_ID}_batch_"

        for part in batch_samplesheets/${EXP_ID}_batch_*; do
            rows_file="\${part}.rows"
            mv "\${part}" "\${rows_file}"
            csv_file="\${part}_samplesheet.csv"
            {
                echo "\$header"
                cat "\${rows_file}"
            } > "\${csv_file}"
            rm -f "\${rows_file}"
        done
    fi

    : > batch_manifest.tsv
    for csv in \$(ls batch_samplesheets/${EXP_ID}_batch_*_samplesheet.csv | sort); do
        batch_id=\$(basename "\$csv" | sed -E 's/^.*_batch_([0-9]{6})_samplesheet\.csv$/batch_\1/')
        printf '%s\t%s\n' "\$batch_id" "\$(realpath "\$csv")" >> batch_manifest.tsv
    done

    echo "Prepared \$(wc -l < batch_manifest.tsv) batch samplesheet(s) for ${EXP_ID}"
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

    container "quay.io/ebigxa/ete3:v2.0"

    containerOptions = "--env BULK_REFERENCES_DIR=${referencePath} --bind ${referencePath}:${referencePath}"

    input:
        path tax_id_file

    output:
        path "*.config"

    script:
    """
    set -euo pipefail

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

    publishDir params.batch_size > 0 ? "${params.outdir}/batches" : "${params.outdir}", mode: 'copy', pattern: "*.rnaseq.done"

    input:
    tuple val(batch_id), path(samplesheet)
    val  EXP_ID
    path star_config
    path params_json, stageAs: "${EXP_ID}_params.json"

    output:
    tuple val(batch_id), path("${batch_id}.rnaseq.done"), emit: done
    tuple val(batch_id), path("${batch_id}.multiqc_report.html"), emit: multiqc_html
    tuple val(batch_id), path("${batch_id}.multiqc_data"), emit: multiqc_data

    script:
    """
    set -euo pipefail

    if [[ ${params.batch_size} -gt 0 ]]; then
        BATCH_OUTDIR="${params.outdir}/batches/${batch_id}"
    else
        BATCH_OUTDIR="${params.outdir}"
    fi
    BATCH_WORKDIR="${nf_workdir}/${EXP_ID}/nested_rnaseq/${batch_id}"
    mkdir -p "\$BATCH_OUTDIR"
    mkdir -p "\$BATCH_WORKDIR"

    echo "Running RNA-seq subworkflow for ${EXP_ID} (${batch_id})"

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
    printf 'STAR_PROFILE_USED\t%s\n' "\${STAR_PROFILE_USED}" > "\$BATCH_OUTDIR/star_profile.log"

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
        --input "${samplesheet}" \\
        --outdir "\$BATCH_OUTDIR" \\
        --contaminant_screening kraken2_bracken \\
        --kraken_db "\${CONTAM_INDEX}" \\
        --without-wave \\
        --save_unaligned \\
        --skip_bbsplit \\
        --skip_pseudo_alignment \\
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
        -work-dir "\$BATCH_WORKDIR" \\
        -with-trace "${params.outdir}/${EXP_ID}_${batch_id}_trace.tsv" \\
        -with-tower \\
        -name "nf_core_rnaseq_${EXP_ID}_${batch_id}" \\
        ${resumeOpt}
    
    # Stage multiqc HTML into task work dir for downstream processes
    MULTIQC_HTML="\$BATCH_OUTDIR/multiqc/star_salmon/multiqc_report.html"
    MULTIQC_DATA_DIR="\$BATCH_OUTDIR/multiqc/star_salmon/multiqc_data"
    if [ ! -f "\$MULTIQC_HTML" ]; then
        echo "ERROR: MultiQC report not found at \$MULTIQC_HTML" >&2
        exit 1
    fi
    if [ ! -d "\$MULTIQC_DATA_DIR" ]; then
        echo "ERROR: MultiQC data directory not found at \$MULTIQC_DATA_DIR" >&2
        exit 1
    fi
    cp "\$MULTIQC_HTML" "${batch_id}.multiqc_report.html"
    cp -R "\$MULTIQC_DATA_DIR" "${batch_id}.multiqc_data"

    # Cleanup nested nf-core work directory only after successful completion.
    if [[ -d "\$BATCH_WORKDIR" ]]; then
        rm -rf "\$BATCH_WORKDIR"
    fi

    # Create done file only if workflow succeeded
    touch "${batch_id}.rnaseq.done"
    
    """
}

process MERGE_BATCH_MULTIQC {

    publishDir "${params.outdir}/multiqc_merged", mode: 'copy', overwrite: true

    conda "bioconda::multiqc=1.27"

    input:
    path multiqc_data_dirs
    val  EXP_ID

    output:
    path "multiqc_report_all_batches.html", emit: merged_html
    path "multiqc_data"
    path "multiqc_merged.done"

    when:
    params.batch_size > 0

    script:
    def mqc_inputs = multiqc_data_dirs.collect { "\"${it}\"" }.join(' ')

    """
    set -euo pipefail

    multiqc \
        -f \
        -o . \
        -n multiqc_report_all_batches.html \
        -c "${workflow.projectDir}/conf/multiqc_star_profile.yaml" \
        ${mqc_inputs}

    touch multiqc_merged.done
    """
}

process MERGE_BATCH_COUNT_MATRICES {

    publishDir "${params.outdir}/merged_counts", mode: 'copy', overwrite: true

    conda "conda-forge::python=3.11 conda-forge::pandas=2.2"

    input:
    path rnaseq_done_files
    val  EXP_ID

    output:
    path "star_salmon/salmon.merged.gene_counts.tsv", emit: gene_counts
    path "star_salmon/salmon.merged.transcript_counts.tsv", emit: transcript_counts
    path "star_salmon/salmon.merged.gene_tpm.tsv", emit: gene_tpm
    path "star_salmon/salmon.merged.transcript_tpm.tsv", emit: transcript_tpm
    path "merged_counts.done", emit: done

    when:
    params.batch_size > 0

    script:
    """
    set -euo pipefail

    mapfile -t batch_ids < <(ls *.rnaseq.done | sed 's/\.rnaseq\.done$//' | sort)
    if [[ \${#batch_ids[@]} -eq 0 ]]; then
        echo "ERROR: No batch completion markers found for ${EXP_ID}" >&2
        exit 1
    fi

    batch_dirs=()
    for batch_id in "\${batch_ids[@]}"; do
        batch_dir="${params.outdir}/batches/\${batch_id}"
        if [[ ! -d "\${batch_dir}" ]]; then
            echo "ERROR: Missing batch output directory: \${batch_dir}" >&2
            exit 1
        fi
        batch_dirs+=("\${batch_dir}")
    done

    python "${workflow.projectDir}/bin/merge_star_salmon_batch_matrices.py" \
        --outdir . \
        "\${batch_dirs[@]}"

    touch merged_counts.done
    """
}

process RECONSTRUCT_BATCH_OUTPUT_LAYOUT {

    publishDir "${params.outdir}", mode: 'copy', overwrite: true, pattern: "layout_reconstructed.done"

    conda "conda-forge::python=3.11 conda-forge::pandas=2.2"

    input:
    path rnaseq_done_files
    val  EXP_ID

    output:
    path "layout_reconstructed.done"

    when:
    params.batch_size > 0

    script:
    """
    set -euo pipefail

    mapfile -t batch_ids < <(ls *.rnaseq.done | sed 's/\.rnaseq\.done$//' | sort)
    if [[ \${#batch_ids[@]} -eq 0 ]]; then
        echo "ERROR: No batch completion markers found for ${EXP_ID}" >&2
        exit 1
    fi

    batch_dirs=()
    for batch_id in "\${batch_ids[@]}"; do
        batch_dir="${params.outdir}/batches/\${batch_id}"
        if [[ ! -d "\${batch_dir}" ]]; then
            echo "ERROR: Missing batch output directory: \${batch_dir}" >&2
            exit 1
        fi
        batch_dirs+=("\${batch_dir}")
    done

    tmp_merge_dir="${params.outdir}/.batch_layout_tmp"
    rm -rf "\${tmp_merge_dir}"

    python "${workflow.projectDir}/bin/merge_nfcore_batch_outputs.py" \
        --outdir "\${tmp_merge_dir}" \
        "\${batch_dirs[@]}"

    # Sync reconstructed content into canonical output directory.
    rsync -a "\${tmp_merge_dir}/" "${params.outdir}/"

    rm -rf "\${tmp_merge_dir}"
    touch layout_reconstructed.done
    """
}

process MATERIALISE_BATCH_FINAL_OUTPUTS {

    publishDir "${params.outdir}/star_salmon", mode: 'copy', overwrite: true, pattern: "salmon.merged.*.tsv"
    publishDir "${params.outdir}/multiqc/star_salmon", mode: 'copy', overwrite: true, pattern: "multiqc_report.html"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, pattern: "final_outputs_materialised.done"

    input:
    path merged_multiqc_html
    path merged_gene_counts
    path merged_tx_counts
    path merged_gene_tpm
    path merged_tx_tpm

    output:
    path "multiqc_report.html"
    path "salmon.merged.gene_counts.tsv"
    path "salmon.merged.transcript_counts.tsv"
    path "salmon.merged.gene_tpm.tsv"
    path "salmon.merged.transcript_tpm.tsv"
    path "final_outputs_materialised.done"

    when:
    params.batch_size > 0

    script:
    """
    set -euo pipefail

    cp "${merged_multiqc_html}" multiqc_report.html
    cp "${merged_gene_counts}" salmon.merged.gene_counts.tsv
    cp "${merged_tx_counts}" salmon.merged.transcript_counts.tsv
    cp "${merged_gene_tpm}" salmon.merged.gene_tpm.tsv
    cp "${merged_tx_tpm}" salmon.merged.transcript_tpm.tsv

    touch final_outputs_materialised.done
    """
}

process SANITISE_MERGED_MULTIQC {

    publishDir "${params.outdir}/multiqc_merged", mode: 'copy', overwrite: true, pattern: "multiqc_report_all_batches*.html"
    publishDir "${params.outdir}/multiqc_merged", mode: 'copy', pattern: "multiqc_sanitisation.done"

    input:
    path merged_multiqc_html

    output:
    path "multiqc_report_all_batches.html", emit: sanitized_html
    path "multiqc_report_all_batches.original.html", emit: original_html
    path "multiqc_sanitisation.done", emit: done

    when:
    params.batch_size > 0

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

    cp "${merged_multiqc_html}" multiqc_report_all_batches.original.html
    sed ${sed_string} multiqc_report_all_batches.original.html > multiqc_report_all_batches.html

    touch multiqc_sanitisation.done
    """
}

process SANITISE_SINGLE_MULTIQC {

    publishDir "${params.outdir}/multiqc/star_salmon", mode: 'copy', overwrite: true, pattern: "multiqc_report.html"
    publishDir "${params.outdir}", mode: 'copy', pattern: "multiqc_sanitisation.done"

    input:
    tuple val(batch_id), path(multiqc_html)

    output:
    path "multiqc_report.html"
    path "multiqc_report.original.html"
    path "multiqc_sanitisation.done"

    when:
    params.batch_size <= 0

    script:
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

    cp "${multiqc_html}" multiqc_report.original.html
    sed ${sed_string} multiqc_report.original.html > multiqc_report.html

    touch multiqc_sanitisation.done
    """
}

process CLEAN_FASTQS {

    publishDir params.outdir, mode: 'copy'

    input:
    path samplesheet
    path rnaseq_done_files

    output:
    path "fastq_cleanup.done"

    script:
    """
    awk -F, 'NR>1 {print \$2; if (\$3 != "") print \$3}' ${samplesheet} \
        | sort -u \
        | xargs -r rm -f

    touch fastq_cleanup.done
    """
}

workflow {
    samplesheet = GET_SAMPLES(params.EXP_ID)

    SET_PARAMS(params.EXP_ID)
    star_config_ch = GET_STAR_PROFILE(SET_PARAMS.out.tax_id)
    params_json_ch = SET_PARAMS.out.params_json

    if (params.batch_size > 0) {
        CREATE_BATCH_SAMPLESHEETS(samplesheet, params.EXP_ID, params.batch_size)
        batch_manifest = CREATE_BATCH_SAMPLESHEETS.out.manifest

        batch_inputs = batch_manifest
            .splitCsv(header: false, sep: '\t')
            .map { row -> tuple(row[0] as String, file(row[1] as String)) }

        RUN_RNASEQ(batch_inputs, params.EXP_ID, star_config_ch, params_json_ch)

        mqc_data_dirs_ch = RUN_RNASEQ.out.multiqc_data.map { it[1] }.collect()
        rnaseq_done_files_ch = RUN_RNASEQ.out.done.map { it[1] }.collect()

        MERGE_BATCH_MULTIQC(mqc_data_dirs_ch, params.EXP_ID)
        MERGE_BATCH_COUNT_MATRICES(rnaseq_done_files_ch, params.EXP_ID)
        RECONSTRUCT_BATCH_OUTPUT_LAYOUT(rnaseq_done_files_ch, params.EXP_ID)
        SANITISE_MERGED_MULTIQC(MERGE_BATCH_MULTIQC.out.merged_html)
        MATERIALISE_BATCH_FINAL_OUTPUTS(
            SANITISE_MERGED_MULTIQC.out.sanitized_html,
            MERGE_BATCH_COUNT_MATRICES.out.gene_counts,
            MERGE_BATCH_COUNT_MATRICES.out.transcript_counts,
            MERGE_BATCH_COUNT_MATRICES.out.gene_tpm,
            MERGE_BATCH_COUNT_MATRICES.out.transcript_tpm
        )
        CLEAN_FASTQS(samplesheet, rnaseq_done_files_ch)
    } else {
        single_batch_inputs = Channel.of(tuple("batch_000001", samplesheet))
        RUN_RNASEQ(single_batch_inputs, params.EXP_ID, star_config_ch, params_json_ch)
        SANITISE_SINGLE_MULTIQC(RUN_RNASEQ.out.multiqc_html)
        rnaseq_done_files_ch = RUN_RNASEQ.out.done.map { it[1] }.collect()
        CLEAN_FASTQS(samplesheet, rnaseq_done_files_ch)
    }
}

workflow.onComplete {
    log.info "Pipeline completed at $workflow.complete."

    if (workflow.success) {
        log.info "Pipeline completed successfully!"

        // def workDir = workflow.workDir
        // if (workDir) {
        //     log.info "Removing work directory: ${workDir}"
        //     workDir.deleteDir()
        // }

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
