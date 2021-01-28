#!/bin/bash

# Renaming files from
# K27ac/BM-k27ac-FBM001_123295_R1.fastq.gz.sort.rmdup.bam
# to
# K27ac/k27ac-FBM001_R1.fastq.gz.sort.rmdup.bam

# Therefore sample ID only unique string. Remove tissue and sequencing ID

for file in *
do
  newname=`echo $file | sed -r "s/[A-Z]{2,3}-//" | sed -r "s/_[0-9]{6}//"`
  mv -i ${file} ${newname}
done

# Had to change file names for my macs2 callpeaks script
