#!/bin/bash

if [ "$#" -ne 3 ]; then
	 echo "Syntax $0 <graphs_dir> <alignments_dir> <calls_dir> <out_dir>"
	 exit 1
fi

GRAPHS=$1
ALIGNMENTS=$2
VARIANTS=$3
OUT_DIR=$4

BASELINE=snp1000g

OPTS="--maxCores 48 --kmer 27 --edge_max 7 --vg_cores 8 --dir_tag"

mkdir $OUT_DIR

#for i in brca1 brca2 sma lrc_kr mhc cenx
for i in brca1 brca2 sma lrc_kir 
do
	 rm -rf blon123 ; ./clusterGraphs.py ./blon123 ${GRAPHS}/*${i}*.vg ${VARIANTS}/${i}/*/*.vg ${OUT_DIR}/${i} ${OPTS} --avg_sample
	 for j in heatmap.pdf heatmap_log.pdf tree.dot  tree.newick  tree.png
	 do
		  cp ${OUT_DIR}/${i}/${j} ${OUT_DIR}/${i}/${i}_avg_${j}
	 done

	 rm -rf blon123 ; ./clusterGraphs.py ./blon123 ${GRAPHS}/*${i}*.vg ${VARIANTS}/${i}/*/*.vg ${OUT_DIR}/${i} ${OPTS}

	 for j in heatmap.pdf heatmap_log.pdf tree.dot  tree.newick  tree.png
	 do
		  cp ${OUT_DIR}/${i}/${j} ${OUT_DIR}/${i}/${i}_${j}
	 done

	 rm -rf blon333 ; ./compareCalls.py ./blon333 ${BASELINE} ${ALIGNMENTS}/${i}/*/*.gam --out_dir ${VARIANTS} --avg_sample --out_sub $i ${OPTS}

	 cp ${VARIANTS}/compare/${i}/*.tsv ${OUT_DIR}/${i}
	 
done
