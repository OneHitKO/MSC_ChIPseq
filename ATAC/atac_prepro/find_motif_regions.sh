#!/bin/bash
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=60G
#SBATCH --time=20:00:00
#SBATCH --output=motif2.stout
#SBATCH --error=motif2.stderr
#SBATCH --export=ALL

source /fast/users/ouk_c/.bashrc
module load gcc
conda init bash
conda activate ng

# Use this script to get the top 25 motifs from unique ATAC peaks w/ or w/o  K27ac 
cd ~/work/msc/ATAC/genrich/peaks/auc200

for value in {1..25}
do
  echo known${value}.motif

  # uniqueBM w/ k27ac
  findMotifsGenome.pl uniqueBM_vs_MSC_K27ac.bed hg19 homer/uniqueBM_vs_MSC_K27ac/ \
    -find homer/uniqueBM_vs_MSC_K27ac/knownResults/known${value}.motif \
    > motif_regions_byK27ac/known${value}.uniqueBM_K27ac.txt
  
  # uniqueBM w/ NO k27ac
  findMotifsGenome.pl uniqueBM_vs_MSC_noK27ac.bed hg19 homer/uniqueBM_vs_MSC_noK27ac/ \
    -find homer/uniqueBM_vs_MSC_noK27ac/knownResults/known${value}.motif \
    > motif_regions_byK27ac/known${value}.uniqueBM_noK27ac.txt

  # uniqueWAT w/ k27ac
  findMotifsGenome.pl uniqueWAT_vs_MSC_K27ac.bed hg19 homer/uniqueWAT_vs_MSC_K27ac/ \
    -find homer/uniqueWAT_vs_MSC_K27ac/knownResults/known${value}.motif \
    > motif_regions_byK27ac/known${value}.uniqueWAT_K27ac.txt

  #findMotifsGenome.pl uniqueWAT_vs_MSC_noK27ac.bed hg19 homer/uniqueWAT_vs_MSC_noK27ac/ \
    -find homer/uniqueWAT_vs_MSC_noK27ac/knownResults/known${value}.motif \
    > motif_regions_byK27ac/known${value}.uniqueWAT_noK27ac.txt

done
