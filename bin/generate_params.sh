#!/usr/bin/env bash

# uses BIOSTUDIES API to fetch species name and utilises genome_reference.conf from bulk-references repo to assign assembly info
# requires genome_reference.conf from bulk-references

usage() { echo """
Usage: 
$0 [ <accession_id> ] 
""" 1>&2; } 

export EXP_ID=$1

BIOSTUDIES_URL="https://www.ebi.ac.uk/biostudies/api/v1/studies/${EXP_ID}"
                        
# Extract species names
species_list=$(curl -fsS "$BIOSTUDIES_URL" \
  | tr -d '\r' \
  | awk '/"name"[[:space:]]*:[[:space:]]*"Organism"/{p=1;next} p&&/"value"/{p=0; sub(/.*"value"[[:space:]]*:[[:space:]]*"/,""); sub(/".*/,""); print}' \
  | sort -u | sed 's/ /_/g' || true)
  
  no=$(printf "%s\n" "${species_list-}" | grep -c . || true)

if [ "$no" -eq 1 ]; then
  printf "%s" "$species_list"   
  export SPECIES=$species_list
else
  >&2 printf "WARN: %s Organism entries for %s\n" "$no" "${EXP_ID}"
  exit 1
fi

echo $SPECIES

genome=$(grep -i ${SPECIES} ../../bulk-references/genome_reference.conf | awk '{print $3}')

if [[ "$genome" == "ensembl" ]]; then
  export RELEASE="$ENSEMBL_RELEASE"
elif [[ "$genome" == "ensemblgenomes" ]]; then
  export RELEASE="$ENSEMBL_GENOME_RELEASE"
else
  echo "$genome is not ensembl or ensemblgenomes"
  exit 1
fi
                                      
export ASSEMBLY=$(grep -i ${SPECIES} ../../bulk-references/genome_reference.conf | awk '{print $7}')

envsubst < "bulk-quantification/params.template.json" > "bulk-quantification/${EXP_ID}_params.json"
