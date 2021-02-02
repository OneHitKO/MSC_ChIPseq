#!/bin/bash
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=60G
#SBATCH --time=20:00:00
#SBATCH --job-name="atac_genrich"
#SBATCH --output=./logs/reproPeaks2_atac_genrich.stdout
#SBATCH --error=./logs/reproPeaks2_atac_genrich.stderr
#SBATCH --export=ALL


# This script takes the ATAC-seq bam files (sorted by readnames in samtools), runs Genrich to call peaks, and ouputs a narrowpeak file
# run from ~/msc/ATAC w/ sbatch 

# Change directory
cd ~/work/msc/ATAC/peaks 

# Output directory
mkdir -p ./reproPeaks

# Setup environment
source ~/.bashrc
conda init bash
conda activate ng
module load gcc

# Paths
bamdir="../bam/sortedQueryname"
hg19="../../deeptools/hg19.blacklist.bed"

## == BM == ##
# get comma separated list of BM bam files
BM=`echo ${bamdir}/BM* | tr ' ' ,`

#Genrich -t $BM \
#  -o reproPeaks/BM_ATAC.narrowPeak \
#  -e chrY -E $hg19 \
#  -j -q 0.05 -a 20 -v
#
## == CH == ##
CH=`echo ${bamdir}/CH* | tr ' ' ,`

Genrich -t $CH \
  -o reproPeaks/CH_ATAC.narrowPeak \
  -e chrY -E $hg19 \
  -j -q 0.05 -a 1 -v 

## == FB == ##
FB=`echo ${bamdir}/FB* | tr ' ' ,`
#
#Genrich -t $FB \
#  -o reproPeaks/FB_ATAC.narrowPeak \
#  -e chrY -E $hg19 \
#  -j -q 0.05 -a 200 -v
#
## == UC == ##
UC=`echo ${bamdir}/UC* | tr ' ' ,`

Genrich -t $UC \
  -o reproPeaks/UC_ATAC.narrowPeak \
  -e chrY -E $hg19 \
  -j -q 0.05 -a 1 -v

## == WAT == ##
WAT=`echo ${bamdir}/WAT* | tr ' ' ,`
#
#Genrich -t $WAT \
#  -o reproPeaks/WAT_ATAC.narrowPeak \
#  -e chrY -E $hg19 \
#  -j -q 0.05 -a 200 -v
#  
