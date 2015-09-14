#!/usr/bin/env bash
# Run after mapping_evaluations.sh or parallelMappingEvaluation.py.
# Makes plots comparing the graphs in each.

set -e

# Grab the input directory to look in
INPUT_DIR=${1}

if [[ ! -d "${INPUT_DIR}" ]]
then
    echo "Specify input directory!"
    exit 1
fi

# Set up the plot parameters
PLOT_PARAMS=(
    --categories "cactus" "camel" "vg" "curoverse" "simons" "snp1000g" "prg" "debruijn-k31" "debruijn-k63"
    "refonly" "trivial" "level1" "level2" "level3"
    --category_labels "Cactus" "Camel" "VG"  "Curoverse" "Simons SNPs" "1000 GSNPs" "PRG" "k=31" "k=63"
    "RefOnly" "Trivial" "Level1" "Level2" "Level3"
    --colors "#5C755E" "#C19A6B" "#000099" "#31184A" "#384DA0" "k" "#353C47" "r" "m" "c" "b" "c" "m" "y"
)

# Where are the stats files
STATS_DIR="${INPUT_DIR}/stats"

# Where do we put the plots?
PLOTS_DIR="${INPUT_DIR}/plots"
mkdir -p "${PLOTS_DIR}"

for REGION in `ls ${STATS_DIR}`
do
    # For every region we ran
    
    # We have intermediate data files for plotting from
    MAPPING_FILE="${PLOTS_DIR}/mapping.${REGION}.tsv"
    MAPPING_PLOT="${PLOTS_DIR}/mapping.${REGION}.png"
    MULTIMAPPING_FILE="${PLOTS_DIR}/multimapping.${REGION}.tsv"
    MULTIMAPPING_PLOT="${PLOTS_DIR}/multimapping.${REGION}.png"
    
    # Make them empty
    :>"${MAPPING_FILE}"
    :>"${MULTIMAPPING_FILE}"
    
    
    echo "Plotting ${REGION^^}..."
    
    for GRAPH_NAME in `ls ${STATS_DIR}/${REGION}`
    do
        # For every graph we ran for it
        
        for STATS_FILE in `ls ${STATS_DIR}/${REGION}/${GRAPH_NAME}`
        do
            # For each sample run, parse its JSON and add a point to the region
            # TSV for the appropriate graph.
            
            # First build the path to the JSON file to look at
            JSON_FILE="${STATS_DIR}/${REGION}/${GRAPH_NAME}/${STATS_FILE}"
            
            # We need to account for the well-mapped/well-multimapped identity thresholds
            
            # First: portion mapped with <=2 mismatches out of 100 expected length
            printf "${GRAPH_NAME}\t" >> "${MAPPING_FILE}"
            # We need the 0 + in case there are no sufficiently good mappings
            cat "${JSON_FILE}" | jq -r '(0 + .primary_mismatches."0" + .primary_mismatches."1" + .primary_mismatches."2") / .total_reads' >> "${MAPPING_FILE}"
            
            # Next: portion NOT multimapped with <=2 mismatches out of 100 expected length
            printf "${GRAPH_NAME}\t" >> "${MULTIMAPPING_FILE}"
            # We need the 0 + in case there are no sufficiently good mappings
            cat "${JSON_FILE}" | jq -r '1 - ((0 + .secondary_mismatches."0" + .secondary_mismatches."1" + .secondary_mismatches."2") / .total_reads)' >> "${MULTIMAPPING_FILE}"
        done
    done
    
    ./boxplot.py "${MAPPING_FILE}" \
        --title "$(printf "Well-mapped (<=2 mismatches)\nreads in ${REGION^^}")" \
        --x_label "Graph" --y_label "Portion Well-mapped" --save "${MAPPING_PLOT}" \
        --x_sideways \
        "${PLOT_PARAMS[@]}" \
        --font_size 20 --dpi 90
        
    ./boxplot.py "${MULTIMAPPING_FILE}" \
        --title "$(printf "Not-well-multimapped\n(>2 mismatches or unmultimapped)\nreads in ${REGION^^}")" \
        --x_label "Graph" --y_label "Portion not-well-multimapped" --save "${MULTIMAPPING_PLOT}" \
        --x_sideways \
        "${PLOT_PARAMS[@]}" \
        --font_size 20 --dpi 90
    
done

