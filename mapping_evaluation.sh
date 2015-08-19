#!/usr/bin/env bash
# mapping_evaluation.sh: do an evaluation of the GA4GH bake-off graphs
# Pipe your server list TSV to the input.
# Format: <region>\t<url>\t<contributor>[\t...] 
# Many copies of this script may be running in parallel.

set -ex

echo "`date`: Starting up on `hostname`..."

# We need to append the server version to the URLs and the input TSV won't have it.
VERSION="v0.6.g"

# Skip the header. TODO: Make sure someone doesn't forget and have their first
# input be silently discarded.
read

mkdir -p alignments
mkdir -p stats
mkdir -p plots
mkdir -p graphs

function run_stats {
    # Call on a basename and a mode (real or sim)
    local BASENAME="${1}"
    local MODE="${2}"

    local STATS_DIR="stats/${REGION}/${MODE}"
    local ALIGNMENT="alignments/${MODE}-${BASENAME}.gam"

    mkdir -p "${STATS_DIR}"

    echo "`date`: Extracting scores..."
    vg view -aj "${ALIGNMENT}" | jq '.score' | grep -v "null" > "${STATS_DIR}/${BASENAME}.scores.tsv"
    
    echo "`date`: Extracting match fractions..."
    # Pull out the reads with non-null scores (i.e. that mapped), total up the
    # length of all the edits that are perfect matches, and spit out the
    # fraction of each read that is a perfect match
    vg view -aj "${ALIGNMENT}" | jq -r 'select(.score != null) | {
        "length": .sequence | length, 
        "matches": ([.path.mapping[].edit[] | select(.to_length == .from_length and .sequence == null) | .to_length] | add)
    } | .matches / .length' > "${STATS_DIR}/${BASENAME}.matches.tsv"
    
    # Work out how many perfect read mappings we have, so we can put it in the plot title
    PERFECT_MAPPINGS=`cat "${STATS_DIR}/${BASENAME}.matches.tsv" | grep '^1$' | wc -l`
    # Keep the number around in a file too.
    echo "${PERFECT_MAPPINGS}" > "${STATS_DIR}/${BASENAME}.perfect.tsv"
    
    echo "`date`: Plotting..."
    local PLOTS_DIR="plots/${MODE}/${REGION}"
    mkdir -p "${PLOTS_DIR}/scores"
    mkdir -p "${PLOTS_DIR}/matches"
    
    # Plot the scores
    ./histogram.py "${STATS_DIR}/${BASENAME}.scores.tsv" \
        --title "${MODE} Mapping Scores for ${BASENAME}" \
        --x_label "Score" \
        --y_label "Read Count" \
        --x_min 0 \
        --bins 50 \
        --save "${PLOTS_DIR}/scores/${BASENAME}.png"
        
    # Plot the perfect match fractions (better since not all reads are the same length)
    ./histogram.py "${STATS_DIR}/${BASENAME}.matches.tsv" \
        --title "${MODE} Mapping Identity (${PERFECT_MAPPINGS} perfect) for ${BASENAME}" \
        --x_label "Match Fraction" \
        --y_label "Read Count" \
        --x_min 0 \
        --bins 50 \
        --save "${PLOTS_DIR}/matches/${BASENAME}.png"
}

while read -r REGION BASE_URL REST
do

    if [[ ${REGION:0:1} == '#' ]]
    then
        # Skip commented-out regions for broken servers
        # See <http://serverfault.com/a/249400> for pound sign weirdness
        continue
    fi
    
    if [[ "${REGION}x" == "x" ]]
    then
        # Skip blank lines
        continue
    fi

    # Since we have namy things with "null" contributor, get the actual
    # submission name from the BASE_URL. It's the part between the last two
    # slashes.
    BASENAME=$(echo "${BASE_URL}" | sed 's/.*\/\(.*\)\/$/\1/')
    echo "`date`: Testing ${BASENAME}"

    # Make a temp dir
    WORK_DIR=$(mktemp -d)

    # Work out the real URL
    URL="${BASE_URL}${VERSION}/"
    
    # Get the graph (in upper case), chop it to remove overly long sequences, and topologically sort and number it.
    echo "`date`: Getting graph..."
    sg2vg "${URL}" -u | vg view -Jv - | vg mod -X 100 - | vg ids -s - > "graphs/${BASENAME}.vg"
    
    # Index it
    echo "`date`: Indexing..."
    rm -rf "graphs/${BASENAME}.vg.index/"
    vg index -s -k10 "graphs/${BASENAME}.vg"
    
    # Do the sim reads
    SIM_FILE="reads/trivial-${REGION^^}.txt"
    
    echo "`date`: Aligning simulated reads..."
    MODE="sim"
    ALIGNMENT="alignments/${MODE}-${BASENAME}.gam"
    
    time vg map -r "${SIM_FILE}" -n 3 "graphs/${BASENAME}.vg" > "${ALIGNMENT}"
    
    run_stats "${BASENAME}" "${MODE}"
    
    # If we have already downloaded reads for the regions, they will be here.
    # Upper-case the region name and add .bam.
    BAM_FILE="reads/${REGION^^}.bam"
    
    echo "`date`: Extracting reads to align..."
    samtools sort -n "${BAM_FILE}" "${WORK_DIR}/reads-by-name"
    bedtools bamtofastq -i "${WORK_DIR}/reads-by-name.bam" -fq "${WORK_DIR}/reads.1.fq" -fq2 "${WORK_DIR}/reads.2.fq" 2>/dev/null
    
    MODE="real"
    ALIGNMENT="alignments/${MODE}-${BASENAME}.gam"
    
    echo "`date`: Aligning real reads..."
    time vg map -f "${WORK_DIR}/reads.1.fq" -f "${WORK_DIR}/reads.2.fq" -n 3 "graphs/${BASENAME}.vg" > "${ALIGNMENT}"
    
    run_stats "${BASENAME}" "${MODE}"
    
    echo "`date`: Surjecting..."
    vg surject -p ref -d "graphs/${BASENAME}.vg.index" -b "alignments/${BASENAME}.gam" > "alignments/${BASENAME}.bam"
    
    # No slashes here so we are protected against variable typos
    rm -rf "${WORK_DIR}"
    
done
