#!/usr/bin/env bash

# Copcy raw fastq files and create a samplesheet.csv file from an input <accession>-configuration.xml or ids.csv.
# The configuration xml file is the definition file for experiments, designed for use in Expression Atlas. 
# The IDs csv file is a single-column list of ENA run IDs, the same input file required by https://nf-co.re/fetchngs/1.12.0/

usage() { echo """
Usage: 
$0 [ -a <accession_id> ] [-x <config.xml>] [-s <samplesheet.csv> ] [-m <era_public_mount_path>] [-c fastq_copy_path ]
or
$0 [ -a <accession_id> ] [-i <ids.csv>] [-s <samplesheet.csv> ] [-m <era_public_mount_path>] [-c fastq_copy_path ]
""" 1>&2; } 

while getopts ":a:x:i:s:m:c:" o; do
    case "${o}" in
        a)
            a=${OPTARG}
            ;;
        x)
            x=${OPTARG}
            ;;
        i)
            i=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        m)
            m=${OPTARG}
            ;;
        c)
            c=${OPTARG}
            ;;
    esac
done
shift $((OPTIND-1))

# Assign and re-assign variables for readability, 

fileIdsType="xml"

if [ -z "${a}" ] || ( [ -z "${x}" ] && [ -z "${i}" ] ) || [ -z "${s}" ] || [ -z "${m}" ]; then
    usage
    exit 1
elif [ -n "${x}" ]; then
    fileIds=$x
else
    fileIds=$i
    fileIdsType="csv"
fi

# Accession has to be an Atlas/BioStudies experiment accession
accession=$a

fileSamples=$s
mountEraPub=$m
copyFastqPath=$c

# Function to derive ENA sub-path from the ENA ID
get_library_subdir() {
    local library=$1
    local forceShortForm=${3:-''} 

    local subDir=${library:0:6}
    local prefix=
    if ! [[ $subDir =~ "ENC" ]] && [[ -z "$forceShortForm" ]] ; then
        local num=${library:3}
        if [ $num -gt 1000000 ]; then

            # ENA pattern is:
            # 
            # - 6-digit codes under e.g. SRR123456
            # - 7-digit codes under e.g. 007/SRR1234567
            # - 8-digit codes under e.g. 078/SRR12345678
            #
            # i.e. we zero-pad to three digits anything after and including the
            # 10th digit. Where we have e.g. '09' we need to strip leading
            # zeros to prevent octal errors with bash.

            digits=$(echo ${library:9} | sed 's/^0*//');
            prefix="$(printf %03d $digits)/"
        fi
    fi
    echo "${subDir}/${prefix}${library}"
}

get_ids_from_input () {
    local fileIds=$1
    local fileIdsType=$2
    if [ $fileIdsType == 'xml' ]; then
        libraries=$( grep "<assay>" "$fileIds" | sed 's/\s*<\/*assay>//g' )
    else
        libraries=$( cat $fileIds )
    fi
    echo "${libraries}"
}

# Main
echo "sample,fastq_1,fastq_2,strandedness" > $fileSamples
for library in $( get_ids_from_input $fileIds $fileIdsType ); do
    librarySubdir=$(get_library_subdir "$library")
    libraryCopyPath="${copyFastqPath}/${accession}"
    mkdir -p $libraryCopyPath

    # Copy subdirectory, without prior knowledge of how many files are inside; should proceed whether or not libraryCopyPath has been created or not 
    aws --no-sign-request --endpoint-url "$FIRE_ENDPOINT" s3 cp "${ERA_PUBLIC_S3_PATH}/${librarySubdir}" "${libraryCopyPath}" --recursive

    libraryFiles=$(find "${libraryCopyPath}" -maxdepth 1 -type f \( -name "${library}*.fastq.gz" -o -name "${library}*.fq.gz" \))

    fileCount=$(echo "${libraryFiles}" | wc -l)
    
    # if [ "$fileCount" -ne 2 ]; then
    #   echo "Error: Expected exactly 2 FASTQ files in ${libraryCopyPath}, but found ${fileCount}."
    #   exit 1
    # fi
    
    # Join the two file paths into a comma-separated list
    libraryFiles=$(echo "${libraryFiles}" | paste -sd "," -)
    
    echo "${library},${libraryFiles},auto" >> $fileSamples
done
