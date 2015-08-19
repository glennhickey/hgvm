#!/usr/bin/env bash

# Get the reads for the bake-off regions
# Requires a samtools with the cram-supporting htslib
# All coordinates are 1-based
# Writes one bam per region, named after the region, in the reads directory.

set -e

# Where does the output go?
OUT_DIR="reads"

# Where are the graphs?
GRAPH_DIR="graphs"

# Point at the correct reference assembly
ALT_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.17_GRCh38.p2/GCA_000001405.17_GRCh38.p2_assembly_structure/all_alt_scaffold_placement.txt"
# We get the alts from there, and infer the regions as the outermost alt endpoints
# TODO: Is that correct?

# Point at the correct sample to get reads from
SAMPLE_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/HG00096/high_cov_alignment/HG00096.alt_bwamem_GRCh38DH.20150424.high_coverage.bam.cram"

# Make the output directory
mkdir -p "${OUT_DIR}"

# Make sure we have the alts file
ALT_FILE=`mktemp`

# Get the alt data from the assembly
echo "Downloading assembly metadata..."
curl "${ALT_URL}" > "${ALT_FILE}" 2>/dev/null

function get_region {
    # For a given region name in any case (like "MHC" or "mhc"), output all the 
    # contig:start-end ranges that we want reads from.
    
    # Grab the region name and upper case it
    local REGION="${1^^}"
    
    case "${REGION}" in
        BRCA1)
            # Not a real alt region. Hardcode the range.
            echo "chr17:43044294-43125482"
            ;;
        BRCA2)
            # Not a real alt region. Hardcode the range.
            echo "chr13:32314861-32399849"
            ;;
        CENX)
            # Not a real alt region. Hardcode this range.
            echo "chrX:58605580-62412542"
            ;;
        *)
            # One of the other alt regions.
            
            # We need this file to hold the relevant alt records
            local REGION_ALTS=`mktemp`
            
            # Get all the alt data and pull out the parent chromosome, alt start
            # and stop on parent, alt name (like GL000254.2), alt start and stop
            # on alt sequence
            cat ${ALT_FILE} | awk "BEGIN {FS=\"\\t\"; OFS=\"\\t\";} {if (\$8 == \"${REGION}\") print \$6,\$12,\$13,\$4,\$10,\$11}" > "${REGION_ALTS}"
            
            # Pull out the chromosome
            CHROMOSOME=`cat "${REGION_ALTS}" | head -n1 | cut -f1`
            # And the start on it
            CHR_START=`cat "${REGION_ALTS}" | cut -f2 | sort -n | head -n1`
            # And the stop on it
            CHR_STOP=`cat "${REGION_ALTS}" | cut -f3 | sort -n | tail -n1`
            
            # Output that
            echo "chr${CHROMOSOME}:${CHR_START}-${CHR_STOP}"
            
            # Rewrite the actual alts names from GL000251.2 into chr6_GL000251v2_alt format and announce them.
            # TODO: This sample doesn't actually have an index covering any alts.
            #cat "${REGION_ALTS}" | cut -f4,5,6 | sed 's/\./v/' | awk "BEGIN {FS=\"\\t\"} {print \"chr${CHROMOSOME}_\" \$1 \"_alt:\" \$2 \"-\" \$3}"
            
            # Clean up the temp file
            rm -f "${REGION_ALTS}"
        ;;
    esac
        
    
}

for REGION in BRCA1 BRCA2 MHC SMA LRC_KIR CENX
do
    echo "Processing ${REGION}..."

    if [[ ! -e "${OUT_DIR}/${REGION}.bam" ]]
    then

        # We need a temp directory to put all our per-range BAMs in
        BAM_DIR=`mktemp -d`

        get_region ${REGION} | 
        {
            # Each range is going to be 0.bam, 1.bam, 2.bam, etc.
            BAM_NUMBER=0
            while read RANGE
            do
                # Download the given range
                echo "Downloading ${RANGE}..."
                samtools view -b -o "${BAM_DIR}/${BAM_NUMBER}.bam" "${SAMPLE_URL}" "${RANGE}"
            
                # Move on to the next range
                BAM_NUMBER=$((BAM_NUMBER + 1))
            done
        }

        if [[ `ls -1 "${BAM_DIR}" | wc -l` == "1" ]]
        then
            echo "Moving ${BAM_NUMBER} bam into ${OUT_DIR}/${REGION}.bam..."
            mv "${BAM_DIR}/0.bam" "${OUT_DIR}/${REGION}.bam"
        else
            # We can only samtools cat with 2 or more files
            echo "Concatenating ${BAM_NUMBER} bams into ${OUT_DIR}/${REGION}.bam..."
            samtools cat -o "${OUT_DIR}/${REGION}.bam" "${BAM_DIR}"/*
        fi
        # Clean up the temporary bams
        rm -rf "${BAM_DIR}"
    else
        echo "${OUT_DIR}/${REGION}.bam already created."
    fi
    
    if [[ "${REGION}" == "CENX" ]]
    then
        # Don't go looking for trivial_cenx because it isn't real.
        # TODO: Use lowest level?
        continue
    fi
    
    # What's the trivial graph?
    TRIVIAL_GRAPH="${GRAPH_DIR}/trivial-${REGION,,}.vg"
    
    # Now get the trivial graph if it doesn't exist
    if [[ ! -e "${TRIVIAL_GRAPH}" ]]
    then
        # Go download it from the trivial URL
        TRIVIAL_URL="http://ec2-54-149-188-244.us-west-2.compute.amazonaws.com/trivial-${REGION,,}/v0.6.g/"
        
        echo "Downloading ${TRIVIAL_GRAPH}..."
        
        # TODO: keep in sync with graph download in mapping_evaluation
        sg2vg "${TRIVIAL_URL}" -u | vg view -Jv - | vg mod -X 100 - | vg ids -s - > "${TRIVIAL_GRAPH}"
    else
        echo "${TRIVIAL_GRAPH} already created..."
    fi
    
    if [[ ! -e "${OUT_DIR}/trivial-${REGION}.txt" ]]
    then
        # Simulate some kmers
        echo "Simulating reads into ${OUT_DIR}/trivial-${REGION}.txt..."
        vg sim -s 1337 -n 10000 -l 250 ${TRIVIAL_GRAPH} > "${OUT_DIR}/trivial-${REGION}.txt"
    else
        echo "${OUT_DIR}/trivial-${REGION}.txt already created"
    fi
done

# What we were doing originally:

#samtools view -b -o BRCA1.bam ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/HG00096/high_cov_alignment/HG00096.alt_bwamem_GRCh38DH.20150424.high_coverage.bam.cram chr17:43044294-43125482

#samtools view -b -o BRCA2.bam ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/HG00096/high_cov_alignment/HG00096.alt_bwamem_GRCh38DH.20150424.high_coverage.bam.cram chr13:32314861-32399849

#samtools view -b -o LRC_KIR.bam ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/HG00096/high_cov_alignment/HG00096.alt_bwamem_GRCh38DH.20150424.high_coverage.bam.cram chr19:54025634-55084318

#samtools view -b -o SMA.bam ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/HG00096/high_cov_alignment/HG00096.alt_bwamem_GRCh38DH.20150424.high_coverage.bam.cram chr5:69216819-71614443

#samtools view -b -o MHC.bam ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/HG00096/high_cov_alignment/HG00096.alt_bwamem_GRCh38DH.20150424.high_coverage.bam.cram chr6:28510120-33480577

#samtools view -b -o CENX.bam ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/HG00096/high_cov_alignment/HG00096.alt_bwamem_GRCh38DH.20150424.high_coverage.bam.cram chrX:58605580-62412542

# Clean up temporary files
rm -f ${ALT_FILE}
