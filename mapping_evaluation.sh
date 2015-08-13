#!/usr/bin/env bash
# mapping_evaluation.sh: do an evaluation of the GA4GH bake-off graphs
# Pipe your server list TSV to the input.
# Format: <region>\t<url>\t<contributor>[\t...] 
set -ex

# We need to append the server version to the URLs and the input TSV won't have it.
VERSION="v0.6.g"

# Skip the header. TODO: Make sure someone doesn't forget and have their first
# input be silently discarded.
read

while read -r REGION BASE_URL CONTRIBUTOR
do

    if [[ ${REGION:0:1} == '#' ]]
    then
        # Skip commented-out regions for broken servers
        # See <http://serverfault.com/a/249400> for pound sign weirdness
        continue
    fi

    # Since we have namy things with "null" contributor, get the actual
    # submission name from the BASE_URL. It's the part between the last two
    # slashes.
    BASENAME=$(echo "${BASE_URL}" | sed 's/.*\/\(.*\)\/$/\1/')
    echo "Testing ${BASENAME}"

    # Make a temp dir
    WORK_DIR=$(mktemp -d)

    # Work out the real URL
    URL="${BASE_URL}${VERSION}/"
    
    # Get the graph
    sg2vg "${URL}" | vg view -Jv - | vg ids -s - > "${WORK_DIR}/${BASENAME}.vg"
    
    # Index it
    vg index -s -k10 "${WORK_DIR}/${BASENAME}.vg"
    
    # If we have already downloaded reads for the regions, they will be here.
    # Upper-case the region name and add .bam.
    BAM_FILE="${REGION^^}.bam"
    
    samtools sort -n "${BAM_FILE}" "${WORK_DIR}/reads-by-name"
    bedtools bamtofastq -i "${WORK_DIR}/reads-by-name.bam" -fq "${WORK_DIR}/reads.1.fq" -fq2 "${WORK_DIR}/reads.2.fq"
    
    vg map -f reads.1.fq -f reads.2.fq "${WORK_DIR}/${BASENAME}.vg" > "${WORK_DIR}/${BASENAME}.gam"
    vg surject -p ref -d "${WORK_DIR}/${BASENAME}.vg.index" -b "${WORK_DIR}/${BASENAME}.gam" > "${WORK_DIR}/${BASENAME}.bam"
    vg view -aj "${WORK_DIR}/${BASENAME}.gam" | jq '.score' | grep -v "null" > "${BASENAME}.scores.tsv"
    
    ./histogram.py "${BASENAME}.scores.tsv" \
        --title "${BASENAME}" \
        --x_label "Score" \
        --y_label "Read Count" \
        --x_min 0 --x_max 210 \
        --bins 20 --line \
        --save "${BASENAME}.png"
        
    # No slashes here so we are protected against variable typos
    rm -rf "${WORK_DIR}"
    
done
