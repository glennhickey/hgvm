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

# We need overall files for mapped and multimapped
OVERALL_MAPPING_FILE="${PLOTS_DIR}/mapping.tsv"
OVERALL_MAPPING_PLOT="${PLOTS_DIR}/mapping.png"
OVERALL_SINGLE_MAPPING_FILE="${PLOTS_DIR}/singlemapping.tsv"
OVERALL_SINGLE_MAPPING_PLOT="${PLOTS_DIR}/singlemapping.png"

for REGION in `ls ${STATS_DIR}`
do
    # For every region we ran
    
    # We have intermediate data files for plotting from
    MAPPING_FILE="${PLOTS_DIR}/mapping.${REGION}.tsv"
    MAPPING_PLOT="${PLOTS_DIR}/mapping.${REGION}.png"
    SINGLE_MAPPING_FILE="${PLOTS_DIR}/singlemapping.${REGION}.tsv"
    SINGLE_MAPPING_PLOT="${PLOTS_DIR}/singlemapping.${REGION}.png"
    RUNTIME_FILE="${PLOTS_DIR}/runtime.${REGION}.tsv"
    RUNTIME_PLOT="${PLOTS_DIR}/runtime.${REGION}.png"
    
    # Make them empty
    :>"${MAPPING_FILE}"
    :>"${SINGLE_MAPPING_FILE}"
    :>"${RUNTIME_FILE}"
    
    
    echo "Plotting ${REGION^^}..."
    
    for GRAPH_NAME in `ls ${STATS_DIR}/${REGION}`
    do
        # For every graph we ran for it
        
        for STATS_FILE in `ls ${STATS_DIR}/${REGION}/${GRAPH_NAME}`
        do
            # For each sample run, parse its JSON and add a point to the region
            # TSV for the appropriate graph.
            
            # We care about single-mapped reads, multi- mapped reads, unmapped
            # reads
            
            # First build the path to the JSON file to look at
            JSON_FILE="${STATS_DIR}/${REGION}/${GRAPH_NAME}/${STATS_FILE}"
            
            # How many reads have any good mapping?
            # We need the 0 + in case there are no sufficiently good mappings
            TOTAL_MAPPED=`cat "${JSON_FILE}" | jq -r '(0 + .primary_mismatches."0" + .primary_mismatches."1" + .primary_mismatches."2")'`
            
            # How many have good secondary mappings? This is a subset of the
            # above.
            TOTAL_MULTIMAPPED=`cat "${JSON_FILE}" | jq -r '(0 + .secondary_mismatches."0" + .secondary_mismatches."1" + .secondary_mismatches."2")'`
            
            # How many reads are there?
            TOTAL_READS=`cat "${JSON_FILE}" | jq -r '.total_reads'`
            
            # Do some math
            PORTION_SINGLE_MAPPED=`echo "(${TOTAL_MAPPED} - ${TOTAL_MULTIMAPPED})/${TOTAL_READS}" | bc -l`
            PORTION_MAPPED=`echo "${TOTAL_MAPPED}/${TOTAL_READS}" | bc -l`
            PORTION_UNMAPPED=`echo "1 - ${TOTAL_MAPPED}/${TOTAL_READS}" | bc -l`
            
            # We need to account for the well-mapped/well-multimapped identity thresholds
            
            # First: portion mapped with <=2 mismatches out of 100 expected length
            printf "${GRAPH_NAME}\t${PORTION_MAPPED}\n" >> "${MAPPING_FILE}"

            # Next: portion single mapped
            printf "${GRAPH_NAME}\t${PORTION_SINGLE_MAPPED}\n" >> "${SINGLE_MAPPING_FILE}"
            
            # Next: runtime in seconds
            printf "${GRAPH_NAME}\t" >> "${RUNTIME_FILE}"
            cat "${JSON_FILE}" | jq -r '.run_time' >> "${RUNTIME_FILE}"
        done
    done
    
    ./boxplot.py "${MAPPING_FILE}" \
        --title "$(printf "Mapped (<=2 mismatches)\nreads in ${REGION^^}")" \
        --x_label "Graph" --y_label "Portion mapped" --save "${MAPPING_PLOT}" \
        --x_sideways \
        "${PLOT_PARAMS[@]}" \
        --font_size 20 --dpi 90
        
    ./boxplot.py "${SINGLE_MAPPING_FILE}" \
        --title "$(printf "Single-mapped (<=2 mismatches)\nreads in ${REGION^^}")" \
        --x_label "Graph" --y_label "Portion single-mapped" --save "${SINGLE_MAPPING_PLOT}" \
        --x_sideways \
        "${PLOT_PARAMS[@]}" \
        --font_size 20 --dpi 90
        
    ./boxplot.py "${RUNTIME_FILE}" \
        --title "$(printf "Aligner runtime\n in ${REGION^^}")" \
        --x_label "Graph" --y_label "Runtime per sample (seconds)" --save "${RUNTIME_PLOT}" \
        --x_sideways \
        "${PLOT_PARAMS[@]}" \
        --font_size 20 --dpi 90
    
done

# Aggregate the overall files
cat "${PLOTS_DIR}"/mapping.*.tsv > "${OVERALL_MAPPING_FILE}"
cat "${PLOTS_DIR}"/singlemapping.*.tsv > "${OVERALL_SINGLE_MAPPING_FILE}"

# Make the overall plots
./boxplot.py "${OVERALL_MAPPING_FILE}" \
    --title "$(printf "Mapped (<=2 mismatches)\nreads")" \
    --x_label "Graph" --y_label "Portion mapped" --save "${OVERALL_MAPPING_PLOT}" \
    --x_sideways \
    "${PLOT_PARAMS[@]}" \
    --font_size 20 --dpi 90

./boxplot.py "${OVERALL_SINGLE_MAPPING_FILE}" \
    --title "$(printf "Single-mapped (<=2 mismatches)\nreads")" \
    --x_label "Graph" --y_label "Portion single-mapped" --save "${OVERALL_SINGLE_MAPPING_PLOT}" \
    --x_sideways \
    "${PLOT_PARAMS[@]}" \
    --font_size 20 --dpi 90

