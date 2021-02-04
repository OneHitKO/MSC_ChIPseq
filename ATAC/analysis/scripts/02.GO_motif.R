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

motif.list[2]


