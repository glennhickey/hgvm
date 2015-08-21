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

    local STATS_DIR="stats/${MODE}/${REGION}"
    local ALIGNMENT="alignments/${MODE}-${BASENAME}.gam"

    mkdir -p "${STATS_DIR}"

    echo "`date`: Extracting best scores..."
    vg view -aj "${ALIGNMENT}" | jq 'select(.is_secondary | not) | .score' | grep -v "null" > "${STATS_DIR}/${BASENAME}.scores.tsv"
    
    echo "`date`: Extracting match fractions..."
    # Pull out the reads with non-null scores (i.e. that mapped), total up the
    # length of all the edits that are perfect matches, and spit out the
    # fraction of each read that is a perfect match in a maximally good alignment for that read
    vg view -aj "${ALIGNMENT}" | jq -r 'select(.score != null and (.is_secondary | not)) | {
        "length": .sequence | length, 
        "matches": ([.path.mapping[].edit[] | select(.to_length == .from_length and .sequence == null) | .to_length] | add)
    } | .matches / .length' > "${STATS_DIR}/${BASENAME}.matches.tsv"
    
    # Work out how many perfect read mappings we have, so we can put it in the plot title
    PERFECT_MAPPINGS=`cat "${STATS_DIR}/${BASENAME}.matches.tsv" | grep '^1$' | wc -l`
    # Keep the number around in a file too.
    echo "${PERFECT_MAPPINGS}" > "${STATS_DIR}/${BASENAME}.perfect.tsv"
    
    echo "`date`: Extracting multimappings..."
    # Since each read can only have 2 mappings, we can count the secondary mappings to get how many reads multimap.
    vg view -aj "${ALIGNMENT}" | jq -r 'select(.score != null and .is_secondary) | {
        "length": .sequence | length, 
        "matches": ([.path.mapping[].edit[] | select(.to_length == .from_length and .sequence == null) | .to_length] | add)
    } | .matches / .length' > "${STATS_DIR}/${BASENAME}.secondary.tsv"
    MULTI_MAPPINGS=`cat "${STATS_DIR}/${BASENAME}.secondary.tsv" | wc -l`
    echo "${MULTI_MAPPINGS}" > "${STATS_DIR}/${BASENAME}.multi.tsv"
    
    echo "`date`: Plotting..."
    local PLOTS_DIR="plots/${MODE}/${REGION}"
    mkdir -p "${PLOTS_DIR}/scores"
    mkdir -p "${PLOTS_DIR}/secondary"
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
        
    # Plot the perfect match fractions for secondary alignments
    ./histogram.py "${STATS_DIR}/${BASENAME}.secondary.tsv" \
        --title "${MODE} Secondary Mapping Identity for ${BASENAME}" \
        --x_label "Match Fraction" \
        --y_label "Read Count" \
        --x_min 0 \
        --bins 50 \
        --save "${PLOTS_DIR}/secondary/${BASENAME}.png"
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
    
    time vg map -r "${SIM_FILE}" -n 3 -M 2 "graphs/${BASENAME}.vg" > "${ALIGNMENT}"
    
    run_stats "${BASENAME}" "${MODE}"
    
    # If we have already downloaded reads for the regions, they will be here.
    # Upper-case the region name and add fastq extensions
    FASTQ1="reads/${REGION^^}.1.fq"
    FASTQ2="reads/${REGION^^}.2.fq"
    
    MODE="real"
    ALIGNMENT="alignments/${MODE}-${BASENAME}.gam"
    
    echo "`date`: Aligning real reads..."
    time vg map -f "${FASTQ1}" -f "${FASTQ1}" -n 3 -M 2 "graphs/${BASENAME}.vg" > "${ALIGNMENT}"
    
    run_stats "${BASENAME}" "${MODE}"
    
    echo "`date`: Surjecting..."
    vg surject -p ref -d "graphs/${BASENAME}.vg.index" -b "alignments/${BASENAME}.gam" > "alignments/${BASENAME}.bam"
    
done
