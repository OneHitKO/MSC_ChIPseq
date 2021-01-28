#!bin/bash

# This script was used to assess insert sizes of ATAC-seq samples

# variable for path to picard.jar
picard=~/work/miniconda/pkgs/picard-2.23.4-0/share/picard-2.23.4-0/picard.jar 

for i in bam/*.bam
do 
  ID=`basename $i .noDUP.noChrM.MAPQ30.bam`
  echo $ID
  
  java -jar $picard CollectInsertSizeMetrics \
    I=${i} \
    O=insertsize/${ID}.txt \
    H=insertsize/${ID}.pdf \
    M=0.5 HISTOGRAM_WIDTH=500
done
