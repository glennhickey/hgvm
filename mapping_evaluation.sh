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

    # Work out the real URL
    URL="${BASE_URL}${VERSION}/"
    
    # Get the graph
    sg2vg "${URL}" | vg view -Jv - | vg ids -s - > "${BASENAME}.vg"
    
    # Index it
    rm -Rf "${BASENAME}.vg.index/"
    vg index -s -k10 "${BASENAME}.vg"
    
    # Get different reads for different regions
    if [[ "${REGION}" == "brca1" ]]
    then
        # BRCA reads from http://web.stanford.edu/~quake/brca/ for now
        # H. Christina Fan et al., Non-invasive Prenatal Measurement of the Fetal Genome, 487 Nature 320 (2012);
        rm -f P028T1-Cell.Free.DNA.BRCA1.bam
        wget http://www.stanford.edu/~quake/brca/P028T1-Cell.Free.DNA.BRCA1.bam
        samtools sort -n P028T1-Cell.Free.DNA.BRCA1.bam reads-by-name
        bedtools bamtofastq -i reads-by-name.bam -fq reads.1.fq -fq2 reads.2.fq
    elif [[ "${REGION}" == "brca2" ]]
    then
        rm -f P028T1-Cell.Free.DNA.BRCA2.bam
        wget http://www.stanford.edu/~quake/brca/P028T1-Cell.Free.DNA.BRCA2.bam
        samtools sort -n P028T1-Cell.Free.DNA.BRCA2.bam reads-by-name
        bedtools bamtofastq -i reads-by-name.bam -fq reads.1.fq -fq2 reads.2.fq
    else 
        echo "No read data available for ${REGION}"
        continue
    fi
    
    vg map -f reads.1.fq -f reads.2.fq "${BASENAME}.vg" > "${BASENAME}.gam"
    vg surject -p ref -d "${BASENAME}.vg.index" -b "${BASENAME}.gam" > "${BASENAME}.bam"
    vg view -aj "${BASENAME}.gam" | jq '.score' | grep -v "null" > "${BASENAME}.scores.txt"
    
    histogram.py "${BASENAME}.scores.txt" \
        --title "${BASENAME}" \
        --x_label "Score" \
        --y_label "Read Count" \
        --x_min 0 --x_max 210 \
        --bins 20 --line \
        --save "${BASENAME}.png"
    
done
