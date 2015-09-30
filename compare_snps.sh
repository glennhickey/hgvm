#!/bin/bash

if [ "$#" -ne 3 ]; then
	 echo "Syntax $0 <graphs_dir> <calls_dir> <out_dir>"
	 exit 1
fi

GRAPHS=$1
VARIANTS=$2
OUT_DIR=$3

OPTS="--maxCores 20 --kmer 27 --edge_max 7"

mkdir -f $OUT_DIR

for i in brca1 brca2 sma lrc_kr mhc
do
	 rm -rf blon123 ; ./clusterGraphs.py ./blon123 ${GRAPHS}/*${i}*.vg ${VARIANTS}/*${i}*.vg ${OUT_DIR}/${i} ${OPTS}
done
