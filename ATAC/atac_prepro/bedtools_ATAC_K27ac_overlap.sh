#!/bin/bash

# Script to get regions w/ ATAC peaks with K27ac to see differential motif and GO term enrichment

# cd to directory containing reproducible unique ATAC peaks
cd ~/work/msc/ATAC/genrich/peaks/auc200

# define location of reproducible k27ac peaks
BM_k27ac=`echo ~/work/msc/CHIP/reproPeaks/BM_K27ac.bed`
WAT_k27ac=`echo ~/work/msc/CHIP/reproPeaks/WAT_K27ac.bed`

# BM unique ATAC peaks w/ and w/o K27ac
bedtools intersect -a uniqueBM_vs_MSC.bed -b $BM_k27ac -v \
  > uniqueBM_vs_MSC_noK27ac.bed

bedtools intersect -a uniqueBM_vs_MSC.bed -b uniqueBM_vs_MSC_noK27ac.bed -v \
  > uniqueBM_vs_MSC_K27ac.bed

# WAT unique ATAC peaks w/ and w/o K27ac
bedtools intersect -a uniqueWAT_vs_MSC.bed -b $WAT_k27ac -v \
  > uniqueWAT_vs_MSC_noK27ac.bed

bedtools intersect -a uniqueWAT_vs_MSC.bed -b uniqueWAT_vs_MSC_noK27ac.bed -v \
  > uniqueWAT_vs_MSC_K27ac.bed
