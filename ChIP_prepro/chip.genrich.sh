#!/bin/bash

# This script will take bam files sorted by samtools -n and identify
# reproducible ChIP seq peaks from biological replicates. 
#
# Genrich also takes into account input files. 
#
# chrY, chrM, and hg19 blacklisted regions are filtered out. 

module load gcc
conda activate ng 

mkdir -p ./reproPeaks

out="./reproPeaks"

### Create comma separated lists for Genrich input!

## define H3K4me3 ##
K4m3BM=`echo resorted/K4m3-BM* | tr ' ' ,`
K4m3CH=`echo resorted/K4m3-CH* | tr ' ' ,`
K4m3FB=`echo resorted/K4m3-FB* | tr ' ' ,`
K4m3UC=`echo resorted/K4m3-UC* | tr ' ' ,`
K4m3WAT=`echo resorted/K4m3-WAT* | tr ' ' ,`

## define H3K27ac ##
K27acBM=`echo resorted/K27ac-BM* | tr ' ' ,`
K27acCH=`echo resorted/K27ac-CH* | tr ' ' ,`
K27acFB=`echo resorted/K27ac-FB* | tr ' ' ,`
K27acUC=`echo resorted/K27ac-UC* | tr ' ' ,`
K27acWAT=`echo resorted/K27ac-WAT* | tr ' ' ,`

## define H3K4me1 ## 
K4m1BM=`echo resorted/K4m1-BM* | tr ' ' ,`
K4m1CH=`echo resorted/K4m1-CH* | tr ' ' ,`
K4m1FB=`echo resorted/K4m1-FB* | tr ' ' ,`
K4m1UC=`echo resorted/K4m1-UC* | tr ' ' ,`
K4m1WAT=`echo resorted/K4m1-WAT* | tr ' ' ,`

## define input ##
inBM=`echo resorted/input-BM* | tr ' ' ,`
inCH=`echo resorted/input-CH* | tr ' ' ,`
inFB=`echo resorted/input-FB* | tr ' ' ,`
inUC=`echo resorted/input-UC* | tr ' ' ,`
inWAT=`echo resorted/input-WAT* | tr ' ' ,`

hg19="../deeptools/hg19.blacklist.bed"
chr="chrY,chrM"

# to do: make more elegant w/ loops

## Call BM peaks
Genrich -t $K4m3BM -c $inBM \
  -o ${out}/BM_K4m3.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 100.0

Genrich -t $K27acBM -c $inBM \
  -o ${out}/BM_K27ac.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 100.0

Genrich -t $K4m1BM -c $inBM \
  -o ${out}/BM_K4m1.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 20.0

## Call CH peaks 
Genrich -t $K4m3CH -c $inCH \
  -o ${out}/CH_K4m3.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 100.0

Genrich -t $K27acCH -c $inCH \
  -o ${out}/CH_K27ac.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 100.0

Genrich -t $K4m1CH -c $inCH \
  -o ${out}/CH_K4m1.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 20.0

## Call FB peaks
Genrich -t $K4m3FB -c $inFB \
  -o ${out}/FB_K4m3.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 100.0

Genrich -t $K27acFB -c $inFB \
  -o ${out}/FB_K27ac.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 100.0

Genrich -t $K4m1FB -c $inFB \
  -o ${out}/FB_K4m1.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 20.0

## Call UC peaks 
Genrich -t $K4m3UC -c $inUC \
  -o ${out}/UC_K4m3.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 100.0

Genrich -t $K27acUC -c $inUC \
  -o ${out}/UC_K27ac.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 100.0

Genrich -t $K4m1UC -c $inUC \
  -o ${out}/UC_K4m1.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 20.0

## Call WAT peaks
Genrich -t $K4m3WAT -c $inWAT \
  -o ${out}/WAT_K4m3.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 100.0

Genrich -t $K27acWAT -c $inWAT \
  -o ${out}/WAT_K27ac.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 100.0

Genrich -t $K4m1WAT -c $inWAT \
  -o ${out}/WAT_K4m1.bed \
  -e $chr -E $hg19 -q 0.05 -v -y -m 30 -a 20.0



















