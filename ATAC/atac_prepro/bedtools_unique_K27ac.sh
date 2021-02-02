#!/bin/bash

# Script to get unique regions to run on homer
# Run script on computing node with ng conda environment

# cd to directory containing reproducible unique ATAC peaks
cd ~/work/msc/ATAC/

# define paths
in="./peaks/reproPeaks"
out="./peaks/uniquePeaks"

# first find strictly unique peaks
bedtools intersect -a ${in}/BM* -b ${in}/UC* ${in}/WAT* ${in}/FB* ${in}/CH* -v
  > ${out}/uniqueBM.narrowPeak 

bedtools intersect -a ${in}/CH* -b ${in}/UC* ${in}/WAT* ${in}/FB* ${in}/BM* -v
  > ${out}/uniqueCH.narrowPeak

bedtools intersect -a ${in}/FB* -b ${in}/BM* ${in}/CH* ${in}/UC* ${in}/WAT* -v
  > ${out}/uniqueFB.narrowPeak 

bedtools intersect -a ${in}/UC* -b ${in}/BM* ${in}/CH* ${in}/FB* ${in}/WAT* -v
  > ${out}/uniqueUC.narrowPeak

bedtools intersect -a ${in}/WAT* -b ${in}/BM* ${in}/CH* ${in}/FB* ${in}/UC* -v
  > ${out}/uniqueWAT.narrowPeak










## define location of reproducible k27ac peaks
#BM_k27ac=`echo ~/work/msc/CHIP/reproPeaks/BM_K27ac.bed`
#WAT_k27ac=`echo ~/work/msc/CHIP/reproPeaks/WAT_K27ac.bed`
#
## BM unique ATAC peaks w/ and w/o K27ac
#bedtools intersect -a uniqueBM_vs_MSC.bed -b $BM_k27ac -v \
#  > uniqueBM_vs_MSC_noK27ac.bed
#
#bedtools intersect -a uniqueBM_vs_MSC.bed -b uniqueBM_vs_MSC_noK27ac.bed -v \
#  > uniqueBM_vs_MSC_K27ac.bed
#
## WAT unique ATAC peaks w/ and w/o K27ac
#bedtools intersect -a uniqueWAT_vs_MSC.bed -b $WAT_k27ac -v \
#  > uniqueWAT_vs_MSC_noK27ac.bed
#
#bedtools intersect -a uniqueWAT_vs_MSC.bed -b uniqueWAT_vs_MSC_noK27ac.bed -v \
#  > uniqueWAT_vs_MSC_K27ac.bed
