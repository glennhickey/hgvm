#!/bin/bash

if [ "$#" -ne 3 ]; then
	 echo "Syntax $0 <graphs_dir> <calls_dir> <out_dir>"
	 exit 1
fi

GRAPHS=$1
ALIGNMENTS=$2
VARIANTS=$3
OUT_DIR=$4

OPTS="--maxCores 24 --kmer 27 --edge_max 7 --vg_cores 2"

mkdir -f $OUT_DIR

for i in brca1 brca2 sma lrc_kr mhc cenx
do
	 rm -rf blan123 ; ./callVariants.py ./blan123 ${ALIGNMENTS}/${i}/*.gam --graph_dir ${GRAPHS} --out_dir ${VARIANTS} ${OPTS}
done
