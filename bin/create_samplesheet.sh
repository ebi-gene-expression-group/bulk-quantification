#!/usr/bin/env bash

# Create a samplesheet.csv file from an input ids.csv

usage() { echo "Usage: $0 [-i <ids.csv>] [-s <samplesheet.csv> ] [-m <era_public_mount_path>]" 1>&2; } 

while getopts ":i:s:m:" o; do
    case "${o}" in
        i)
            i=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        m)
            m=${OPTARG}
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${i}" ] || [ -z "${s}" ] || [ -z "${m}" ]; then
    usage
    exit 1
fi

# Re-assign variables for readability
fileIds=$i
fileSamples=$s
mountEraPub=$m

# Function to derive ENA sub-path from the ENA ID
get_library_path() {
    local library=$1
    local rootDir=$2
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
    echo "${mountEraPub}/${subDir}/${prefix}${library}"
}


# Main
echo "sample,fastq_1,fastq_2,strandedness" > $fileSamples
while IFS= read -r library; do
    libraryPath=$(get_library_path  "$library" "$mountEraPub")
    libraryFiles=$(find "${libraryPath}" -maxdepth 1 -type f | paste -sd "," - ) 
    echo "${libraryFiles},auto" >> $fileSamples
done < "$fileIds"