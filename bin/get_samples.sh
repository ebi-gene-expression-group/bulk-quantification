#!/usr/bin/env bash

# Copy raw fastq files and create a samplesheet.csv file.
# Either provide a file listing run IDs <ids.csv>, or state that input is from atlas.
# The IDs csv file is a single-column list of ENA run IDs, the same input file required by https://nf-co.re/fetchngs/1.12.0/

usage() { echo """
Usage: 
$0 [ -a <accession_id> ] [ -x atlas ] [ -s <samplesheet.csv> ] [ -e <endpoint_url> ] [ -m <era_public_mount_path or era_public_s3_path>] [ -c fastq_copy_path ]
or
$0 [ -a <accession_id> ] [ -x <ids.csv> ] [ -s <samplesheet.csv> ] [ -e <endpoint_url> ] [-m <era_public_mount_path or era_public_s3_path>] [-c fastq_copy_path ]
""" 1>&2; } 

while getopts ":a:x:s:e:m:c:" o; do
    case "${o}" in
        a)
            a=${OPTARG}
            ;;
        x)
            x=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        e)
            e=${OPTARG}
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

if [ -z "${a}" ] || [ -z "${x}" ] || [ -z "${s}" ] || [ -z "${e}" ] || [ -z "${m}" ]; then
    usage
    exit 1
fi

# Accession has to be an Atlas/BioStudies experiment accession
accession=$a
fileIds=$x
fileSamples=$s
endpointUrl=$e
eraPubPath=$m
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
    local accession=$1
    local fileIds=$2
    if [ $fileIds == 'atlas' ]; then
        fileConfigXml=$(ls ${ATLAS_PROD}/analysis/*/rna-seq/experiments/${accession}/${accession}-configuration.xml)
        libraries=$( grep "</assay>" "${fileConfigXml}" | sed -n 's/.*<assay[^>]*>\(.*\)<\/assay>.*/\1/p' | sort -u )
    else
        libraries=$( cat $fileIds )
    fi
    echo "${libraries}"
}

# Main
for library in $( get_ids_from_input $accession $fileIds ); do
    echo "library id to be downloaded $library"
    librarySubdir=$(get_library_subdir "$library")
    echo "library subdir $librarySubdir"
    libraryCopyPath="${copyFastqPath}/${accession}"
    echo "library CopyPath $libraryCopyPath"
    mkdir -p $libraryCopyPath

    # Copy subdirectory, without prior knowledge of how many files are inside; should proceed whether or not libraryCopyPath has been created or not 
    aws --no-sign-request --endpoint-url "${endpointUrl}" s3 cp "${eraPubPath}/${librarySubdir}" "${libraryCopyPath}" --recursive

    # First look for paired-end files
    pairedFiles=$(find "${libraryCopyPath}" -maxdepth 1 -type f \
      -name "${library}_[12].f*q.gz")
    
    if [[ -n "$pairedFiles" ]]; then
      libraryFiles="$pairedFiles"
    else
      libraryFiles=$(find "${libraryCopyPath}" -maxdepth 1 -type f \
        -name "${library}.f*q.gz")
    fi

    fileCount=$(echo "${libraryFiles}" | wc -l)
    echo "fileCount $fileCount"
    if [ ! -s "$fileSamples" ]; then
        if [ "$fileCount" -eq 1 || "$fileCount" -eq 2 ]; then
            echo "sample,fastq_1,fastq_2,strandedness" > "$fileSamples"
        else
            echo "Error: Expected exactly 1 or 2 FASTQ files in ${libraryCopyPath}, but found ${fileCount}."
            exit 1
        fi
    fi
    
    # Join the two file paths into a comma-separated list
    libraryFiles=$(echo "${libraryFiles}" | paste -sd "," -)

    [[ "$libraryFiles" != *,* ]] && libraryFiles="${libraryFiles},"
    
    echo "${library},${libraryFiles},auto" >> $fileSamples
done
