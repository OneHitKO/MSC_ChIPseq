#!/bin/bash
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=60G
#SBATCH --time=20:00:00
#SBATCH --job-name="get regions"
#SBATCH --output=logs/knownRegions.stdout
#SBATCH --error=logs/knownRegions.stderr
#SBATCH --export=ALL

# This script retrieves the peak IDs of de novo Homer motifs identified from
# unique ATAC seq peaks
# Use as sbatch from msc/ATAC folder (so all logs in one place)

# config
source /fast/users/ouk_c/.bashrc
module load gcc
conda init bash
conda activate ng

# set paths
cd ~/work/msc/ATAC/
bed="./peaks/uniquePeaks"

# loop through all homer results
for results in ./results_homer/*/
do
 # create new folder with peaks
 mkdir -p ${results}/homerRegions/
 
 # get cell type name of homer results folder
 cell=`echo $results | cut -d'/' -f 3`
  
 # find regions in top 5 de novo motifs
 for value in {1..5}
 do
   findMotifsGenome.pl ${bed}/${cell}.narrowPeak hg19 $results \
     -find ${results}/homerResults/motif${value}.motif \
     > ${results}/homerRegions/${cell}_motif${value}.txt
 done

done

