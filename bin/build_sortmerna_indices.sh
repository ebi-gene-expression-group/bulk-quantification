#!/bin/bash
set -euo pipefail

# Configuration
REFERENCE_BASE="${BULK_REFERENCES_DIR}"
GENOME_VERSION="GRCh38"
SILVA_VERSION="138.2"

# Paths
GENOME_DIR="${REFERENCE_BASE}/homo_sapiens/Ensembl/${GENOME_VERSION}/Sequence/WholeGenomeFasta/Homo_sapiens.${GENOME_VERSION}.dna.toplevel.fa"
GENE_DIR="${REFERENCE_BASE}/homo_sapiens/Ensembl/${GENOME_VERSION}/Annotation/Genes/Homo_sapiens.${GENOME_VERSION}.114.gtf"
RIBO_DIR="${REFERENCE_BASE}/contamination/ribo"
INDEX_DIR="${REFERENCE_BASE}/contamination/sortmerna/index"

# Create a temporary working directory
WORK_DIR="${NF_WORKDIR}/ribo_index/"
trap "rm -rf ${WORK_DIR}" EXIT

echo "=== Building SortMeRNA Index ==="
echo "Genome: ${GENOME_VERSION}"
echo "SILVA: ${SILVA_VERSION}"
echo "Output: ${INDEX_DIR}"
echo ""

# Create minimal samplesheet for index building only
cat > ${WORK_DIR}/dummy_samplesheet.csv << EOF
sample,fastq_1,fastq_2,strandedness
dummy,${WORK_DIR}/dummy_R1.fq.gz,${WORK_DIR}/dummy_R2.fq.gz,auto
EOF

# Create minimal dummy fastq files (just for passing validation)
echo "@read1" | gzip > ${WORK_DIR}/dummy_R1.fq.gz
echo "ACGT" | gzip >> ${WORK_DIR}/dummy_R1.fq.gz
echo "+" | gzip >> ${WORK_DIR}/dummy_R1.fq.gz
echo "IIII" | gzip >> ${WORK_DIR}/dummy_R1.fq.gz

cp ${WORK_DIR}/dummy_R1.fq.gz ${WORK_DIR}/dummy_R2.fq.gz

# Run nf-core/rnaseq with only index building
nextflow run ${workflow.projectDir}/subworkflows/rnaseq/main.nf \
    --input ${WORK_DIR}/dummy_samplesheet.csv \
    --outdir ${WORK_DIR}/output \
    --fasta ${GENOME_DIR} \
    --gtf ${GENE_DIR} \
    --remove_ribo_rna \
    --ribo_removal_tool sortmerna \
    --ribo_database_manifest ${RIBO_DIR}/silva.manifest.tsv \
    --save_reference \
    --skip_alignment \
    --skip_pseudo_alignment \
    --skip_fastqc \
    --skip_trimming \
    -profile singularity \
    -work-dir ${WORK_DIR}/work \
    -resume

# Copy index to final location
echo ""
echo "Copying index to final location..."
mkdir -p ${INDEX_DIR}
cp -rv ${WORK_DIR}/output/genome/index/sortmerna/* ${INDEX_DIR}/

# Verify index
if ls ${INDEX_DIR}/*.stats >/dev/null 2>&1; then
    echo "✓ SortMeRNA index successfully created!"
    echo "  Location: ${INDEX_DIR}"
    echo "  Files:"
    ls -lh ${INDEX_DIR}
    
    # Create a metadata file
    cat > ${INDEX_DIR}/index.metadata << EOF
creation_date: $(date -Iseconds)
genome_version: ${GENOME_VERSION}
silva_version: ${SILVA_VERSION}
manifest: ${RIBO_DIR}/silva.manifest.tsv
nf_core_rnaseq_version: 3.14.0
EOF
    
else
    echo "✗ Error: Index creation failed!"
    exit 1
fi

echo ""
echo "=== Index Building Complete ==="
