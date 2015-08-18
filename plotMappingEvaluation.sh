#!/usr/bin/env bash
# Run after mapping_evaluations.sh or parallelMappingEvaluation.py.
# Makes plots comparing the graphs in each.

set -e

COLLATED_DIR="stats/collated"
mkdir -p "${COLLATED_DIR}"

PLOTS_DIR="plots/collated"
mkdir -p "${PLOTS_DIR}"

for MODE in real sim
do
    mkdir -p "${COLLATED_DIR}/${MODE}"
    mkdir -p "${PLOTS_DIR}/${MODE}"
    for REGION in brca1 brca2 mhc sma lrc_kir cenx
    do
        echo "Processing ${REGION} mode ${MODE}..."
    
        COLLATED_FILE="${COLLATED_DIR}/${MODE}/${REGION}.tsv"
        rm -f "${COLLATED_FILE}"
        rm -f "${PLOTS_DIR}/${MODE}/${REGION}.png"
        
        # Make sure we get an empty list of files if none exist.
        shopt -s nullglob
        
        for FILE in stats/${REGION}/${MODE}/*.perfect.tsv
        do
            # Put the category name with no newline
            # Make sure to strip out the region name (no matter which side its dash is on)
            printf "${FILE}\t" | sed 's/.*\/\(.*\)\.perfect\.tsv/\1/' | \
                sed "s/-${REGION}//g" | sed "s/${REGION}-//g" >> "${COLLATED_FILE}"
        
            # Put the number of perfect reads in each category
            cat "${FILE}" >> "${COLLATED_FILE}"
        done
        
        if [[ "${MODE}" == "sim" ]]
        then
            # For simulations we want to set a max
            MAX="10000"
        else
            # For real reads we want the read count of the input bam
            MAX=`samtools flagstat reads/${REGION^^}.bam | head -n1 | sed 's/\([0-9]*\).*/\1/'`
        fi
        
        if [[ -e "${COLLATED_FILE}" ]]
        then
        
            # Now make the actual plot
            # Here we specify the order and colors for the different schemes
            # Only ones with data are used
            ./barchart.py "${COLLATED_FILE}" --title "Perfectly mapped ${MODE} reads in ${REGION^^}" \
                --x_label "Graph" --y_label "Read Count" --save "${PLOTS_DIR}/${MODE}/${REGION}.png" \
                --x_sideways --min 0 --max "${MAX}" \
                --categories "cactus" "camel" "debruijn-k31" "debruijn-k63" "refonly" "trivial" \
                "level1" "level2" "level3" \
                --colors "g" "y" "r" "m" "c" "b" "c" "m" "y"
        fi
    done
done
