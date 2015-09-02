#!/usr/bin/env bash
# Run after mapping_evaluations.sh or parallelMappingEvaluation.py.
# Makes plots comparing the graphs in each.

set -e

# What's the minimum identity to consider a read well-mapped or well-
# multimapped?
WELL_MAPPED_THRESHOLD=0.98

COLLATED_DIR="stats/collated"
mkdir -p "${COLLATED_DIR}"

PLOTS_DIR="plots/collated"
mkdir -p "${PLOTS_DIR}"

for MODE in real sim
do
    mkdir -p "${COLLATED_DIR}/${MODE}"
    mkdir -p "${PLOTS_DIR}/${MODE}"
    
    rm -f "${COLLATED_DIR}/${MODE}"/all-*.tsv
    rm -f "${PLOTS_DIR}/${MODE}"/all-*.png
    
    for REGION in brca1 brca2 mhc sma lrc_kir cenx
    do
        echo "Processing ${REGION} mode ${MODE}..."
    
        DATA_FILE="${COLLATED_DIR}/${MODE}/${REGION}-primary-perfect.tsv"
        rm -f "${DATA_FILE=}"
        PLOT_FILE="${PLOTS_DIR}/${MODE}/perfect-${REGION}.png"
        rm -f "${PLOT_FILE}"
        
        # Make sure we get an empty list of files if none exist.
        shopt -s nullglob
        
        for FILE in stats/${MODE}/${REGION}/*.perfect.tsv
        do
            # Put the category name with no newline
            # Make sure to strip out the region name (no matter which side its dash is on)
            printf "${FILE}\t" | sed 's/.*\/\(.*\)\.perfect\.tsv/\1/' | \
                sed "s/-${REGION}//g" | sed "s/${REGION}-//g" >> "${DATA_FILE}"
        
            # Put the number of perfect reads in each category
            cat "${FILE}" >> "${DATA_FILE}"
            
            # TODO: Normalize for use in the all file
        done
        
        if [[ -e "${DATA_FILE}" ]]
        then
        
            # Now make the actual plot
            # Here we specify the order and colors for the different schemes
            # Only ones with data are used
            ./barchart.py "${DATA_FILE}" --title "$(printf "Perfectly-mapped ${MODE}\nreads in ${REGION^^}")" \
                --x_label "Graph" --y_label "Read Count" --save "${PLOT_FILE}" \
                --categories "cactus" "camel" "vg" "curoverse" "simons" "snp1000g" "prg" "debruijn-k31" "debruijn-k63"\
                "refonly" "trivial" "level1" "level2" "level3" \
                --category_labels "Cactus" "Camel" "VG"  "Curoverse" "Simons SNPs" "1000G SNPs" "PRG" "k=31" "k=63"\
                "RefOnly" "Trivial" "Level1" "Level2" "Level3" \
                --colors "#5C755E" "#C19A6B" "#000099" "#31184A" "#384DA0" "k" "#353C47" "r" "m" "c" "b" "c" "m" "y" \
                --x_sideways \
                --font_size 20 --dpi 90
        fi
        
        # How many reads were tried?
        if [[ "${MODE}" == "sim" ]]
        then
            # For simulations we want to set a max
            MAX="10000"
        else
            # For real reads we want the read count of the input bam
            MAX=`cat reads/${REGION^^}.*.fq | grep "^@" | wc -l`
        fi
        
        DATA_FILE="${COLLATED_DIR}/${MODE}/${REGION}-primary-good.tsv"
        rm -f "${DATA_FILE}"
        PLOT_FILE="${PLOTS_DIR}/${MODE}/good-${REGION}.png"
        rm -f "${PLOT_FILE}"
        
        # Make sure we get an empty list of files if none exist.
        shopt -s nullglob
        
        for FILE in stats/${MODE}/${REGION}/*.matches.tsv
        do
            # What graph is this?
            BASENAME=`echo "${FILE}" | sed 's/.*\/\(.*\)\.matches\.tsv/\1/'`
            
            # Work out how many primary alignments are good enough
            OVER_COUNT=`cat "${FILE}" | awk "{if (\\\$1 >= \${WELL_MAPPED_THRESHOLD}) print \\\$0}" | wc -l`
            # And how many reads exist
            TOTAL_COUNT="${MAX}"
            
            # Put the category name with no newline in the file to plot
            # Make sure to strip out the region name (no matter which side its dash is on)
            printf "${FILE}\t" | sed 's/.*\/\(.*\)\.matches\.tsv/\1/' | \
                sed "s/-${REGION}//g" | sed "s/${REGION}-//g" >> "${DATA_FILE}"
            
            # Put the fraction of reads where the primary alignments are sufficiently good.
            echo "${OVER_COUNT} / ${TOTAL_COUNT}" | bc -l >> "${DATA_FILE}"
            
            # Put a point in the all file, assigned to this graph.
            printf "${FILE}\t" | sed 's/.*\/\(.*\)\.matches\.tsv/\1/' | \
                sed "s/-${REGION}//g" | sed "s/${REGION}-//g" >> "${COLLATED_DIR}/${MODE}/all-good.tsv"
            echo "${OVER_COUNT} / ${TOTAL_COUNT}" | bc -l >> "${COLLATED_DIR}/${MODE}/all-good.tsv"
            
        done
        
        if [[ -e "${DATA_FILE}" ]]
        then
        
            # Now make the actual plot
            # Here we specify the order and colors for the different schemes
            # Only ones with data are used
            ./barchart.py "${DATA_FILE}" \
                --title "$(printf "Well-mapped (>=${WELL_MAPPED_THRESHOLD}) ${MODE}\nreads in ${REGION^^}")" \
                --x_label "Graph" --y_label "Portion Well-mapped" --save "${PLOT_FILE}" \
                --categories "cactus" "camel" "vg" "curoverse" "simons" "snp1000g" "prg" "debruijn-k31" "debruijn-k63"\
                "refonly" "trivial" "level1" "level2" "level3" \
                --category_labels "Cactus" "Camel" "VG"  "Curoverse" "Simons SNPs" "1000G SNPs" "PRG" "k=31" "k=63"\
                "RefOnly" "Trivial" "Level1" "Level2" "Level3" \
                --colors "#5C755E" "#C19A6B" "#000099" "#31184A" "#384DA0" "k" "#353C47" "r" "m" "c" "b" "c" "m" "y" \
                --x_sideways \
                --font_size 20 --dpi 90
        fi
        
        # Now what portion of secondary reads aren't good?
        DATA_FILE="${COLLATED_DIR}/${MODE}/${REGION}-secondary-bad.tsv"
        rm -f "${DATA_FILE}"
        PLOT_FILE="${PLOTS_DIR}/${MODE}/secondary-${REGION}.png"
        rm -f "${PLOT_FILE}"
        
        for FILE in stats/${MODE}/${REGION}/*.secondary.tsv
        do
            # What graph is this?
            BASENAME=`echo "${FILE}" | sed 's/.*\/\(.*\)\.secondary\.tsv/\1/'`
            
            # Work out how many secondary alignments are bad enough
            UNDER_COUNT=`cat "${FILE}" | awk "{if (\\\$1 < \${WELL_MAPPED_THRESHOLD}) print \\\$0}" | wc -l`
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
                sed "s/-${REGION}//g" | sed "s/${REGION}-//g" >> "${DATA_FILE}"
            
            # Put the fraction of reads where the secondary alignments are sufficiently bad.
            echo "${UNDER_COUNT} / ${TOTAL_COUNT}" | bc -l >> "${DATA_FILE}"
            
            # Put a point in the all file, assigned to this graph.
            printf "${FILE}\t" | sed 's/.*\/\(.*\)\.secondary\.tsv/\1/' | \
                sed "s/-${REGION}//g" | sed "s/${REGION}-//g" >> "${COLLATED_DIR}/${MODE}/all-secondary.tsv"
            echo "${UNDER_COUNT} / ${TOTAL_COUNT}" | bc -l >> "${COLLATED_DIR}/${MODE}/all-secondary.tsv"
            
        done
        
        if [[ -e "${DATA_FILE}" ]]
        then
        
            # Now make the actual plot
            # Here we specify the order and colors for the different schemes
            # Only ones with data are used
            ./barchart.py "${DATA_FILE}" \
                --title "$(printf "Non-well-multimapped (<${WELL_MAPPED_THRESHOLD}) ${MODE}\nreads in ${REGION^^}")" \
                --x_label "Graph" --y_label "$(printf "Portion Not\nMultimapped")" --save "${PLOT_FILE}" \
                --categories "cactus" "camel" "vg" "curoverse" "simons" "snp1000g" "prg" "debruijn-k31" "debruijn-k63"\
                "refonly" "trivial" "level1" "level2" "level3" \
                --category_labels "Cactus" "Camel" "VG" "Curoverse" "Simons SNPs" "1000G SNPs" "PRG" "k=31" "k=63"\
                "RefOnly" "Trivial" "Level1" "Level2" "Level3" \
                --colors "#5C755E" "#C19A6B" "#000099" "#31184A" "#384DA0" "k" "#353C47" "r" "m" "c" "b" "c" "m" "y" \
                --x_sideways \
                --font_size 20 --dpi 90
        fi
    done
    
    if [[ -e "${COLLATED_DIR}/${MODE}/all-good.tsv" ]]
    then
        # Plot a dot plot of the portion of reads with good primary alignments
        ./scatter.py "${COLLATED_DIR}/${MODE}/all-good.tsv" \
            --dotplot \
            --title "$(printf "Well-mapped (>=${WELL_MAPPED_THRESHOLD}) ${MODE}\nreads across regions")" \
            --x_label "Graph" --y_label "Portion Well-Mapped" --save "${PLOTS_DIR}/${MODE}/all-good.png" \
            --categories "cactus" "camel" "vg" "curoverse" "simons" "snp1000g" "prg" "debruijn-k31" "debruijn-k63"\
            "refonly" "trivial" \
            --colors "#5C755E" "#C19A6B" "#000099" "#31184A" "#384DA0" "k" "#353C47" "r" "m" "c" "b" \
            --min_y 0 --max_y 1.05 \
            --font_size 10 --dpi 90
    fi
    
    if [[ -e "${COLLATED_DIR}/${MODE}/all-secondary.tsv" ]]
    then
        # Plot a dot plot of the portion of reads with no good secondary alignments
        ./scatter.py "${COLLATED_DIR}/${MODE}/all-secondary.tsv" \
            --dotplot \
            --title "$(printf "Non-well-multimapped (<${WELL_MAPPED_THRESHOLD}) ${MODE}\nreads across regions")" \
            --x_label "Graph" --y_label "Portion Not Multimapped" --save "${PLOTS_DIR}/${MODE}/all-secondary.png" \
            --categories "cactus" "camel" "vg" "curoverse" "simons" "snp1000g" "prg" "debruijn-k31" "debruijn-k63"\
            "refonly" "trivial" \
            --colors "#5C755E" "#C19A6B" "#000099" "#31184A" "#384DA0" "k" "#353C47" "r" "m" "c" "b" \
            --min_y 0 --max_y 1.05 \
            --font_size 10 --dpi 90
    fi
done
