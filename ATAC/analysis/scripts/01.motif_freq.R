#!/usr/bin/env Rscript

##--- Analyzes and plots the frequency of known motifs identified in unique ATAC peaks ---## 
# run from ATAC/analysis/scripts directory under conda r environment

# load libraries
library("tidyverse")

##-- create list of "known motif" homer results of unique peaks --##
dir = "~/work/msc/ATAC/results_homer/"

files = list.files(path = dir, pattern = "knownResults.txt", recursive = T,
                  full.names = T)

motif.list = map(files, read_tsv)

# use stringr regexp to extract "uniqueBM", etc 
names(motif.list) = str_extract(files, pattern = "unique[:upper:]{2,3}")

#-- data reformating --# 
# create new column to label which cells the results came from 
motif.list = mapply(cbind, motif.list, "ATAC" = names(motif.list), 
                    SIMPLIFY = F) 

# rename cols, join tibbles
motif.list = map(motif.list, ~ dplyr::rename(., motif = 1,
                                             p_val = 3,
                                             log_p_val = 4,
                                             q_val = 5,
                                             num_targets = 6,
                                             perc_targets = 7,
                                             num_bg = 8,
                                             perc_bg = 9))

motif.tbl = do.call("rbind", motif.list)

# coerce perc cols to double to calculate ratio over bg
motif.tbl$perc_targets = as.double(str_replace(motif.tbl$perc_targets, 
                                               pattern = "%", replacement = ""))

motif.tbl$perc_bg = as.double(str_replace(motif.tbl$perc_bg,
                                          pattern = "%", replacement = ""))

tail(motif.tbl)










