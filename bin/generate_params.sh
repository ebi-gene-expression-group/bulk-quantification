#!/usr/bin/env bash

# uses BIOSTUDIES API to fetch species name and utilises genome_reference.conf from bulk-references repo to assign assembly info
# requires genome_reference.conf from bulk-references

usage() { echo """
Usage: 
$0 [ <accession_id> ] 
""" 1>&2; } 

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export EXP_ID=$1

# Function to fetch species names from BioStudies API
fetch_species_names() {
  local exp_id="$1"
  local BIOSTUDIES_URL="https://www.ebi.ac.uk/biostudies/api/v1/studies/${exp_id}"
  local species_list no
  local exp_type=$(echo "${exp_id#E-}" | sed 's/-.*//')

    species_list="$(
      awk -F'\t' '
        NR==1 {
          for (i=1; i<=NF; i++) if ($i == "Characteristics [organism]" || $i == "Characteristics[organism]") col=i
          if (!col) { print "ERROR: column Characteristics [organism] not found" > "/dev/stderr"; exit 1 }
          next
        }
        { sub(/\r$/, "", $col); print $col }
      ' "$AE2_PRODUCTION/${exp_type}/${exp_id}/${exp_id}.sdrf.txt" \
      | sort -u \
      | sed 's/ /_/g'
    )"
  

  no="$(printf '%s\n' "${species_list}" | grep -c . || true)"

  if [[ "$no" -eq 1 ]]; then
    local SP="$species_list"
    echo "$SP" >&2
  else
    >&2 printf "WARN: %s Organism entries for %s\n" "$no" "$exp_id"
    return 1
  fi
}


export SPECIES=$(fetch_species_names "${EXP_ID}")
export SPECIES_lower="$(tr '[:upper:]' '[:lower:]' <<<"$SPECIES")"

genome=$(grep -i ${SPECIES} $SCRIPT_DIR/../../bulk-references/genome_reference.conf | awk '{print $3}')
tax_id=$(grep -i ${SPECIES} $SCRIPT_DIR/../../bulk-references/genome_reference.conf | awk '{print $2}')
RELEASE=""
if [[ "$genome" == "ensembl" ]]; then
  RELEASE="$ENSEMBL_RELEASE"
elif [[ "$genome" == "ensemblgenomes" ]]; then
  RELEASE="$ENSEMBL_GENOME_RELEASE"
else
  echo "$genome is not ensembl or ensemblgenomes" >&2
  exit 1
fi

export RELEASE
                                      
export ASSEMBLY=$(grep -i ${SPECIES} $SCRIPT_DIR/../../bulk-references/genome_reference.conf | awk '{print $7}')

envsubst < "$SCRIPT_DIR/../params.template.json" > "${EXP_ID}_params.json"

echo $tax_id
