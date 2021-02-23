#!/bin/bash
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=60G
#SBATCH --time=20:00:00
<<<<<<< HEAD
#SBATCH --job-name="atac_msc"
#SBATCH --output=./logs/repro_atac_genrich.stdout
#SBATCH --error=./logs/repro_atac_genrich.stderr
=======
#SBATCH --job-name="atac_genrich"
#SBATCH --output=./logs/reproPeaks2_atac_genrich.stdout
#SBATCH --error=./logs/reproPeaks2_atac_genrich.stderr
>>>>>>> prepro
#SBATCH --export=ALL


# This script takes the ATAC-seq bam files (sorted by readnames in samtools), runs Genrich to call peaks, and ouputs a narrowpeak file
<<<<<<< HEAD
# The ideal number of peaks should be at least > 50,000 peaks 

# Change directory
cd ~/work/msc/ATAC/peaks
=======
# run from ~/msc/ATAC w/ sbatch 

# Change directory
cd ~/work/msc/ATAC/peaks 
>>>>>>> prepro

# Output directory
mkdir -p ./reproPeaks

<<<<<<< HEAD
# Setup environment 
source ~/.bashrc

=======
# Setup environment
source ~/.bashrc
>>>>>>> prepro
conda init bash
conda activate ng
module load gcc

<<<<<<< HEAD
# Paths 
bamdir="../bam/sortedQueryname"
hg19="../../deeptools/hg19.blacklist.bed"

# Call peaks with Genrich (combining biological replicates to call reproducible peaks)
# Script is updated to increase peaks called by changing AUC is cell-specific manner
=======
# Paths
bamdir="../bam/sortedQueryname"
hg19="../../deeptools/hg19.blacklist.bed"
>>>>>>> prepro

## == BM == ##
# get comma separated list of BM bam files
BM=`echo ${bamdir}/BM* | tr ' ' ,`
<<<<<<< HEAD

Genrich -t $BM \
  -o reproPeaks/BM_ATAC.narrowPeak \
  -e chrY -E $hg19 \
  -j -q 0.05 -a 25 -v
=======
>>>>>>> prepro

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
<<<<<<< HEAD
  -j -q 0.05 -a 10 -v 

## == FB == ##
#FB=`echo ${bamdir}/FB* | tr ' ' ,`
=======
  -j -q 0.05 -a 1 -v 

## == FB == ##
FB=`echo ${bamdir}/FB* | tr ' ' ,`
>>>>>>> prepro
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
<<<<<<< HEAD
  -j -q 0.05 -a 10 -v

## == WAT == ##
#WAT=`echo ${bamdir}/WAT* | tr ' ' ,`
=======
  -j -q 0.05 -a 1 -v

## == WAT == ##
WAT=`echo ${bamdir}/WAT* | tr ' ' ,`
>>>>>>> prepro
#
#Genrich -t $WAT \
#  -o reproPeaks/WAT_ATAC.narrowPeak \
#  -e chrY -E $hg19 \
#  -j -q 0.05 -a 200 -v
#  
