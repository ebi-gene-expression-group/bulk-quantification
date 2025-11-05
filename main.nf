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

// Define output directory based on EXP_ID and ATLAS_PROD
params.outdir = "${atlasProd}/test_analysis/baseline/rna-seq/experiments/${params.EXP_ID}"

workflow {
    samplesheet = create_samplesheet(params.EXP_ID)
    species_ch = GET_SPECIES(params.EXP_ID)
    run_rnaseq(samplesheet, params.EXP_ID, species_ch)
}


// -------------------- PROCESSES -------------------- //

process create_samplesheet {
    publishDir "${params.outdir}/samplesheet", mode: 'copy'

    input:
    val EXP_ID

    output:
    path "${EXP_ID}_samplesheet.csv"

    script:
    """
    echo "Creating samplesheet for ${EXP_ID}"

    CONFIG_FILE="${params.outdir}/${EXP_ID}-configuration.xml"
    SAMPLESHEET="${EXP_ID}_samplesheet.csv"

    bash "${projectDir}/bin/create_samplesheet.sh" -x "\${CONFIG_FILE}" -s "\${SAMPLESHEET}" -m "${era_public_mount_path}"
    """
}


process GET_SPECIES {
  input:
    val EXP_ID

  output:
    stdout emit: species

  script:
  """
  set -euo pipefail

  URL="https://www.ebi.ac.uk/biostudies/api/v1/studies/${EXP_ID}"

  # Extract species names
  species_list=\$(curl -fsS "\$URL" \
    | tr -d '\\r' \
    | awk '/"name"[[:space:]]*:[[:space:]]*"Organism"/{p=1;next} p&&/"value"/{p=0; sub(/.*"value"[[:space:]]*:[[:space:]]*"/,""); sub(/".*/,""); print}' \
    | sort -u | sed 's/ /_/g' || true)

  no=\$(printf "%s\\n" "\${species_list-}" | grep -c . || true)

  if [ "\$no" -eq 1 ]; then
    printf "%s" "\$species_list"   # stdout → val species
  else
    >&2 printf "WARN: %s Organism entries for %s\\n" "\$no" "$EXP_ID"
    exit 1
  fi
  """
}


process run_rnaseq {
    publishDir "${params.outdir}/rnaseq", mode: 'copy'

    input:
    path samplesheet
    val  EXP_ID
    val SPECIES

    output:
    path "${EXP_ID}.rnaseq.done"

    script:
    """
    set -euo pipefail

    export EXP_ID="${EXP_ID}"
    export SAMPLESHEET="${samplesheet}"
    export OUTDIR="${params.outdir}/rnaseq"

    export SPECIES="${SPECIES}"

    echo "Running RNA-seq subworkflow for \${EXP_ID} (species=\${SPECIES})"

    # Render params file from template (must reference \$EXP_ID, \$SAMPLESHEET, \$OUTDIR, \$SPECIES)
    envsubst < "${projectDir}/params.template.json" > "\${EXP_ID}_params.json"
    echo "Rendered params:"
    cat "\${EXP_ID}_params.json"

    nextflow run ${projectDir}/subworkflows/rnaseq/main.nf \\
        -params-file "\${EXP_ID}_params.json" \\
        -c "${projectDir}/conf/rnaseq.config" \\
        -profile singularity \\
        --without-wave

    # Success flag for downstream logic / idempotency
    touch "\${EXP_ID}.rnaseq.done"
    """
}
