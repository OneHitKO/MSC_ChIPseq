#!/bin/bash

# This script takes the ATAC-seq bam files (sorted by readnames in samtools), runs Genrich to call peaks,
# and ouputs a narrowpeak file
#
# Change directory
cd ~/work/msc/ATAC/genrich

# Output directory
#mkdir -p ./peaks

# Setup environment 
conda activate ng

# Call peaks with Genrich (combining biological replicates to call reproducible peaks)

# to do: make more elegant w/ loops

## == BM == ##
# get comma separated list of BM bam files
BM=`echo resort/BM* | tr ' ' ,`

Genrich -t $BM \
  -o peaks/BM_ATAC_AUC50.narrowPeak \
  -e chrY -E ~/work/msc/deeptools/hg19.blacklist.bed \
  -j -q 0.05 -a 200 -v

## == CH == ##
CH=`echo resort/CH* | tr ' ' ,`

Genrich -t $CH \
  -o peaks/CH_ATAC_AUC50.narrowPeak \
  -e chrY -E ~/work/msc/deeptools/hg19.blacklist.bed \
  -j -q 0.05 -a 200 -v 

## == FB == ##
FB=`echo resort/FB* | tr ' ' ,`

Genrich -t $FB \
  -o peaks/FB_ATAC_AUC50.narrowPeak \
  -e chrY -E ~/work/msc/deeptools/hg19.blacklist.bed \
  -j -q 0.05 -a 200 -v

## == UC == ##
UC=`echo resort/UC* | tr ' ' ,`

Genrich -t $UC \
  -o peaks/UC_ATAC_AUC50.narrowPeak \
  -e chrY -E ~/work/msc/deeptools/hg19.blacklist.bed \
  -j -q 0.05 -a 200 -v

## == WAT == ##
WAT=`echo resort/WAT* | tr ' ' ,`

Genrich -t $WAT \
  -o peaks/WAT_ATAC_AUC50.narrowPeak \
  -e chrY -E ~/work/msc/deeptools/hg19.blacklist.bed \
  -j -q 0.05 -a 200 -v
  
