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

# create fold over bg column 
motif.tbl$fold = motif.tbl$perc_targets/motif.tbl$perc_bg

# reformat motif column to motif, family (study not necessary)
motif.tbl = separate(motif.tbl, motif, into = c("Motif","Family"),
                     sep = "\\(|\\)")


# practice plot
p1 = filter(motif.tbl, fold >= 2, q_val < 0.01) %>%
ggplot(., aes(fold, -(log_p_val)))+
       geom_point(aes(color = Family, size = perc_targets))+
       theme_bw()+
       facet_wrap(vars(ATAC), scales = "free_y")+
       theme(axis.text.y = element_text(size = 8),
             panel.grid = element_blank(),
             legend.box = "vertical")+
       scale_y_continuous(labels = function(x) format(x,scientific = T))

ggsave(plot = p1, filename = "../plots/motif_freq.png", height = 7, width = 7)










