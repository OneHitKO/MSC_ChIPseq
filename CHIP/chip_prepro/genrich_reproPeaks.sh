#!/bin/bash

# This script will take bam files sorted by samtools -n and identify
# reproducible ChIP seq peaks from biological replicates. 
#
# Genrich also takes into account input files. 
#
# chrY, chrM, and hg19 blacklisted regions are filtered out. 

cd ~/work/msc/CHIP/

#module load gcc
#conda init bash
#conda activate ng 

# input dir with resorted bam files
in="./resorted/"

# create output dir
mkdir -p ./peaks/reproPeaks
out="./peaks/reproPeaks"

# other param for Genrich
hg19="../deeptools/hg19.blacklist.bed"
chr="chrY,chrM"

# create tissue array
tissue=("BM" "CH" "FB" "UC" "WAT")


# loop through each tissue
for i in "${tissue[@]}"
do
  echo $i

  # fine resorted bam files matching tissue pattern
  k27ac=`echo $(find ${in}/K27ac*.bam -type f -print | grep "\-${i}\-")`
  k4m3=`echo $(find ${in}/K4m3*.bam -type f -print | grep "\-${i}\-")`
  k4m1=`echo $(find ${in}/K4m1*.bam -type f -print | grep "\-${i}\-")`

  
  # perform peak calling separately (might need to adjust param)
  # chip-seq do NOT have PE reads! use -y argument!
  echo "$k27ac"
  Genrich -t "${k27ac}" -c "${k27ac//K27ac/input}" \
    -o ${out}/${i}_k27ac.narrowPeak \
    -e $chr -E $hg19 -y -m 30 \
    -q 0.05 -a 100.0

  echo "$k4m3"
  Genrich -t "${k4m3}" -c "${k4m3//K4m3/input}" \
    -o ${out}/${i}_k4m3.narrowPeak \
    -e $chr -E $hg19 -y -m 30 \
    -q 0.05 -a 100.0

  echo "$k4m1"
  Genrich -t "${k4m1}" -c "${k4m1//K4m1/input}" \
    -o ${out}/${i}_k4m1.narrowPeak \
    -e $chr -E $hg19 -y -m 30 \
    -q 0.1 -a 50.0 -g 400

done

