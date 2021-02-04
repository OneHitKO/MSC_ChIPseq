#!/usr/bin/env Rscript

###--- Script will perform GO term analysis on motif family-specific regions ---###

library(tidyverse)
library(rtracklayer)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(org.Hs.eg.db)
#library(ChIPseeker)
#library(clusterProfiler)

#-- import unique ATAC peaks --# 
dir = "../../peaks/uniquePeaks/"

bed.files = list.files(path=dir, pattern=".narrowPeak$", full.names=T)

# create list of granges of unique peaks
peaks.list = purrr::map(bed.files, ~rtracklayer::import(., format="narrowPeak"))

names(peaks.list) = str_extract(bed.files, pattern = "unique[:upper:]{2,3}")

#-- import motif files --# 
dir = "../../results_homer/"

# create nested list containing info on peaks enriched in motif familys
motif.list = vector("list",5)

names(motif.list) = names(peaks.list)

for (i in seq_along(motif.list)){
  
  # get vector of file paths to regions for each motif family
  motif.files = list.files(path = paste0(dir, names(motif.list)[i],"/homerRegions/"), 
                                         full.names = T, pattern = ".txt$")
  
  # import files and create nested list 
  motif.list[[i]] = purrr::map(motif.files, read_tsv)

  # create vector to get motif names to name nexted list
  motifnames = vector("character", length(motif.list[[i]]))

  for (j in seq_along(motif.list[[i]])){
    
    # add name of motif by getting a matrix of split strings
    motifnames[j] = str_split_fixed(motif.list[[i]][[j]]$`Motif Name`,
                                    pattern = ":|/", n = 2)[1,2]
  }
           
  # rename nested list
  names(motif.list[[i]]) = motifnames
}

#-- subset granges which contain specific motifs --#

GR.motifs = vector("list", length = 5)
names(GR.motifs) = names(motif.list)

for (i in seq_along(GR.motifs)){
  
  # initialize the nested list
  GR.motifs[[i]] = vector("list", length(motif.list[[i]]))
  names(GR.motifs[[i]]) = names(motif.list[[i]])

  # subset granges, probably a better way to do this
  for (j in seq_along(motif.list[[i]])){
    GR.motifs[[i]][[j]] = peaks.list[[i]][mcols(peaks.list[[i]])$name %in% motif.list[[i]][[j]]$PositionID]
   } 
}

#-- subset motif-specific, cell-specifc granges with k27ac --# 

# first import reproducible narrow peaks (called by genrich)
dir = "../../../CHIP/peaks/reproPeaks/" 
 
k27ac.files = list.files(path=dir, pattern="_K27ac.bed$", recursive=T,
                       full.names=T)

k27ac.list = purrr::map(k27ac.files, ~ rtracklayer::import(.,format="narrowPeak"))

names(k27ac.list)=str_extract(k27ac.files, pattern ="[:upper:]{2,3}" )

names(k27ac.list)



#-- combine +k27ac and -k27ac lists --#




#-- annotate --# 




#-- create function --# 
