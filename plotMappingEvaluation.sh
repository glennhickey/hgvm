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
    
        COLLATED_PRIMARY_FILE="${COLLATED_DIR}/${MODE}/${REGION}-primary-perfect.tsv"
        rm -f "${COLLATED_PRIMARY_FILE}"
        rm -f "${PLOTS_DIR}/${MODE}/primary-${REGION}.png"
        
        # Make sure we get an empty list of files if none exist.
        shopt -s nullglob
        
        for FILE in stats/${MODE}/${REGION}/*.perfect.tsv
        do
            # Put the category name with no newline
            # Make sure to strip out the region name (no matter which side its dash is on)
            printf "${FILE}\t" | sed 's/.*\/\(.*\)\.perfect\.tsv/\1/' | \
                sed "s/-${REGION}//g" | sed "s/${REGION}-//g" >> "${COLLATED_PRIMARY_FILE}"
        
            # Put the number of perfect reads in each category
            cat "${FILE}" >> "${COLLATED_PRIMARY_FILE}"
        done
        
        # Now what portion of secondary reads aren't good?
        COLLATED_SECONDARY_FILE="${COLLATED_DIR}/${MODE}/${REGION}-secondary-bad.tsv"
        rm -f "${COLLATED_SECONDARY_FILE}"
        rm -f "${PLOTS_DIR}/${MODE}/secondary-${REGION}.png"
        
        for FILE in stats/${MODE}/${REGION}/*.secondary.tsv
        do
            # What graph is this?
            BASENAME=`echo "${FILE}" | sed 's/.*\/\(.*\)\.secondary\.tsv/\1/'`
            
            # Work out how many secondary alignments are bad enough
            UNDER_COUNT=`cat "${FILE}" | awk '{if ($1 < 0.98) print $0}' | wc -l`
            # And how many secondary ones
            MULTI_COUNT=`cat "${FILE}" |  wc -l`
            # And how many primary ones exist
            TOTAL_COUNT=`cat stats/${MODE}/${REGION}/${BASENAME}.scores.tsv | wc -l`
            
            # Add the count of primary ones with no secondary mappings to the
            # count of sufficiently bad secondary mappings.
            UNDER_COUNT=$((UNDER_COUNT + TOTAL_COUNT - MULTI_COUNT))
            
            # Put the category name with no newline in the file to plot
            # Make sure to strip out the region name (no matter which side its dash is on)
            printf "${FILE}\t" | sed 's/.*\/\(.*\)\.secondary\.tsv/\1/' | \
                sed "s/-${REGION}//g" | sed "s/${REGION}-//g" >> "${COLLATED_SECONDARY_FILE}"
            
            # Put the fraction of reads where the secondary alignments are sufficiently bad.
            echo "${UNDER_COUNT} / ${TOTAL_COUNT}" | bc -l >> "${COLLATED_SECONDARY_FILE}"
            
            # TODO: nonexistent isn't sufficiently bad
        done
        
        if [[ "${MODE}" == "sim" ]]
        then
            # For simulations we want to set a max
            MAX="10000"
        else
            # For real reads we want the read count of the input bam
            MAX=`samtools flagstat reads/${REGION^^}.bam | head -n1 | sed 's/\([0-9]*\).*/\1/'`
        fi
        
        if [[ -e "${COLLATED_PRIMARY_FILE}" ]]
        then
        
            # Now make the actual plot
            # Here we specify the order and colors for the different schemes
            # Only ones with data are used
            ./barchart.py "${COLLATED_PRIMARY_FILE}" --title "$(printf "Perfectly mapped ${MODE}\nreads in ${REGION^^}")" \
                --x_label "Graph" --y_label "Read Count" --save "${PLOTS_DIR}/${MODE}/primary-${REGION}.png" \
                --categories "cactus" "camel" "vg" "curoverse" "simons" "snp1000g" "prg" "debruijn-k31" "debruijn-k63"\
                "refonly" "trivial" "level1" "level2" "level3" \
                --category_labels "Cactus" "Camel" "VG"  "Curoverse" "Simons SNPs" "1000G SNPs" "PRG" "k=31" "k=63"\
                "RefOnly" "Trivial" "Level1" "Level2" "Level3" \
                --colors "#5C755E" "#C19A6B" "#000099" "#31184A" "#384DA0" "k" "#353C47" "r" "m" "c" "b" "c" "m" "y" \
                --x_sideways \
                --font_size 20 --dpi 90
        fi
        
        if [[ -e "${COLLATED_SECONDARY_FILE}" ]]
        then
        
            # Now make the actual plot
            # Here we specify the order and colors for the different schemes
            # Only ones with data are used
            ./barchart.py "${COLLATED_SECONDARY_FILE}" --title "$(printf "Non-multimapped ${MODE}\nreads in ${REGION^^}")" \
                --x_label "Graph" --y_label "$(printf "Portion Not\nMultimapped")" --save "${PLOTS_DIR}/${MODE}/secondary-${REGION}.png" \
                --categories "cactus" "camel" "vg" "curoverse" "simons" "snp1000g" "prg" "debruijn-k31" "debruijn-k63"\
                "refonly" "trivial" "level1" "level2" "level3" \
                --category_labels "Cactus" "Camel" "VG" "Curoverse" "Simons SNPs" "1000G SNPs" "PRG" "k=31" "k=63"\
                "RefOnly" "Trivial" "Level1" "Level2" "Level3" \
                --colors "#5C755E" "#C19A6B" "#000099" "#31184A" "#384DA0" "k" "#353C47" "r" "m" "c" "b" "c" "m" "y" \
                --x_sideways \
                --font_size 20 --dpi 90
        fi
    done
done
