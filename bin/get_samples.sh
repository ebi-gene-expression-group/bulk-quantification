#!/usr/bin/env bash

# Copy raw fastq files and create a samplesheet.csv file.
# Either provide a file listing run IDs <ids.csv>, or state that input is from atlas.
# The IDs csv file is a single-column list of ENA run IDs, the same input file required by https://nf-co.re/fetchngs/1.12.0/

# WARNING: This version currently only works for `-x atlas`

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

        :)  echo "ERROR: Option -$OPTARG requires an argument." >&2
            usage
            ;;

        \?) echo "ERROR: Invalid option: -$OPTARG" >&2
            usage
            ;;
    
    esac
done

shift $((OPTIND-1))

# Assign and re-assign variables for readability, 
if [ -z "${a}" ] || [ -z "${x}" ] || [ -z "${s}" ] || [ -z "${e}" ] || [ -z "${m}" ]; then
    echo "ERROR: Missing argument(s)." >&2
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

validate_assay_id() {
    local assay_id=$1

    if [[ ! "$assay_id" =~ ^([DES]RR[0-9]+|ENC[A-Za-z0-9_.-]+)$ ]]; then
        echo "ERROR: Unsafe or unsupported assay ID: ${assay_id}" >&2
        echo "       Expected an ERR/SRR/DRR run accession or an ENCODE identifier starting with ENC." >&2
        exit 1
    fi
}

get_ids_from_input () {
    local accession=$1
    local fileIds=$2
    local fileConfigXml
    local libraries

    if [[ $fileIds == atlas ]]; then
        shopt -s nullglob
        local matches=( "${ATLAS_PROD}/analysis/"*/rna-seq/experiments/"${accession}/${accession}-configuration.xml" )

        if (( ${#matches[@]} != 1 )); then
            echo "ERROR: expected exactly one config XML for ${accession}, found ${#matches[@]}" >&2
            exit 1
        fi

        fileConfigXml=${matches[0]}

        libraries=$(
            grep '</assay>' "$fileConfigXml" |
            sed -n 's/.*<assay[^>]*>\(.*\)<\/assay>.*/\1/p' |
            sort -u
        )
    else
        libraries=$(<"$fileIds")
    fi

    printf '%s\n' "$libraries"
}

# Main
while IFS= read -r library; do
    validate_assay_id "$library"

    echo "Library ID to be downloaded: ${library}"

    librarySubdir=$(get_library_subdir "$library")
    echo "Library subdir: ${librarySubdir}"

    libraryCopyPath="${copyFastqPath}/${accession}"
    echo "Library copy path: ${libraryCopyPath}"
    mkdir -p "$libraryCopyPath"

    # List files to download, ordered by filename
    expectedFiles=$( aws --no-sign-request --endpoint-url "${endpointUrl}" s3 ls "${eraPubPath}/${librarySubdir}/" | awk '{  print $4 }' | sort )
    fileCount=$( echo "${expectedFiles}" | wc -w )

    # Exit if there are more than 2 files in ENA
    if (( fileCount != 1 && fileCount != 2 )); then
        echo "ERROR: Expected exactly 1 or 2 FASTQ files in ENA, but found ${fileCount}."
        echo "       Files: ${expectedFiles}"
        exit 1
    fi

    # Download files
    aws --no-sign-request --endpoint-url "${endpointUrl}" s3 sync "${eraPubPath}/${librarySubdir}" "${libraryCopyPath}" --no-progress

    # Check if all expected files were copied
    # To include in the future: md5sum validation (but this will need checking the ENA db)
    dlExit=false
    dlFiles=""
    for libFile in $expectedFiles; do
        dlFiles+="${libraryCopyPath}/${libFile} "
        if [ ! -s "${libraryCopyPath}/${libFile}" ]; then
            echo "ERROR: File not downloaded properly: ${libFile}"
            dlExit=true
        fi
    done

    # Exit with error if not all files were downloaded successfully
    if $dlExit; then
        echo "Exiting because not all files were downloaded properly for ${library}."
        exit 1
    fi

    # Write entry in the samplesheet
    if [ ! -s "$fileSamples" ]; then
        echo "sample,fastq_1,fastq_2,strandedness" > "$fileSamples"
    fi
        
    # Join the file path/s into a comma-separated list
    libraryFiles=$(printf '%s' "$dlFiles" | sed 's/\s\+/,/g')
    if [[ $fileCount -eq 1 ]]; then
        libraryFiles="${libraryFiles},"
    fi
    echo "${library},${libraryFiles},auto" >> "$fileSamples"

done < <(get_ids_from_input "$accession" "$fileIds")
