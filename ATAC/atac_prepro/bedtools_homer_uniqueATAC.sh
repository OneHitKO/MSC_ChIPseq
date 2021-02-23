#!/bin/bash

# Script to identify enriched motifs on unique ATAC-seq Peaks (they have not
# been stratified by H3K27ac)
# Run script on computing node with ng conda environment

# cd to directory containing reproducible unique ATAC peaks
cd ~/work/msc/ATAC/

# define paths
in="./peaks/reproPeaks"
out="./peaks/uniquePeaks"

# first find strictly unique peaks
bedtools intersect -a ${in}/BM* -b ${in}/UC* ${in}/WAT* ${in}/FB* ${in}/CH* -v > ${out}/uniqueBM.narrowPeak 

bedtools intersect -a ${in}/CH* -b ${in}/UC* ${in}/WAT* ${in}/FB* ${in}/BM* -v > ${out}/uniqueCH.narrowPeak

bedtools intersect -a ${in}/FB* -b ${in}/BM* ${in}/CH* ${in}/UC* ${in}/WAT* -v > ${out}/uniqueFB.narrowPeak 

bedtools intersect -a ${in}/UC* -b ${in}/BM* ${in}/CH* ${in}/FB* ${in}/WAT* -v > ${out}/uniqueUC.narrowPeak

bedtools intersect -a ${in}/WAT* -b ${in}/BM* ${in}/CH* ${in}/FB* ${in}/UC* -v > ${out}/uniqueWAT.narrowPeak


# motif enrichment analysis on unique ATAC seq peaks
for i in ${out}/*.narrowPeak
do
  id=`basename $i .narrowPeak`

  findMotifsGenome.pl $i hg19 ./results_homer/${id}/ \
    -size 200 -S 10 -mask
done
