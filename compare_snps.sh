#!/bin/bash

for i in brca1 brca2 sma lrc_kr brc
do
	 rm -rf blon ; ./compareGraphs.py ./blon graphs/snp1000g-${i}.vg variants/kmer/real/*${i}*.vg compare/${i}_snp1000_comp.tsv
done
