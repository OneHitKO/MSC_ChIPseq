#!/bin/bash

# This script takes the bam files from histone ChIP seq, runs macs2 callpeak and outputs bed files for downstream analysis.
# Usage sh macs2.sh

# Change directory
cd ~/work/msc/bam

# Output directories
mkdir -p ./macs2output

# Set up environment
conda activate macs2

# Find ID names, call peaks
for i in K27ac/*.bam
do

  # Extract ID from .bam files
  sample=`echo $i | sed -r 's/[kK][0-9]{1,2}[a-m0-9]{1,2}-//'`
  ID=`basename $sample _R1.fastq.gz.sort.rmdup.bam`

  # Call K27ac narrowPeaks
  echo "calling K27ac peaks for $ID"

  macs2 callpeak -t K27ac/k27ac-${ID}_R1.fastq.gz.sort.rmdup.bam \
  -c input/input-${ID}_R1.fastq.gz.sort.rmdup.bam \
  -f BAM \
  -g hs \
  -n ${ID}_K27ac \
  --outdir ./macs2output/K27ac

  # Call K4m3 narrowPeaks
  echo "calling K4m3 peaks for $ID"

  macs2 callpeak -t K4m3/K4m3-${ID}_R1.fastq.gz.sort.rmdup.bam \
  -c input/input-${ID}_R1.fastq.gz.sort.rmdup.bam \
  -f BAM \
  -g hs \
  -n ${ID}_K4m3 \
  --outdir ./macs2output/K4m3

  # Call broad marks for K4m1
  echo "calling K4m1 peaks for $ID"

  macs2 callpeak -t K4m1/K4m1-${ID}_R1.fastq.gz.sort.rmdup.bam \
  -c input/input-${ID}_R1.fastq.gz.sort.rmdup.bam \
  --broad \
  -f BAM \
  -g hs \
  -n ${ID}_K4m1 \
  --outdir ./macs2output/K4m1

  echo "complete!"
done
