#!/usr/bin/env Rscript

###--- Script will perform GO term analysis on unique peaks---###

library(tidyverse)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
library(clusterProfiler)

#-- import unique ATAC peaks --# 
dir = "../uniquePeaks/"

bed.files = list.files(path=dir, pattern=".narrowPeak$", full.names=T)

# create list of granges of unique peaks
peaks.list = purrr::map(bed.files, ~rtracklayer::import(., format="narrowPeak"))

names(peaks.list) = str_extract(bed.files, pattern = "unique[:upper:]{2,3}")

#-- add cell info in mcols --#
# define function to use in map
add_mcol = function(gr_obj,col_name,value){
  mcols(gr_obj)[[toString(col_name)]] = value
  return(gr_obj)
}

peaks.list = purrr::map2(peaks.list, names(peaks.list), ~ add_mcol(.x, "Tissue", .y))

#-- overlap w/ k27ac info --#
k27ac.list = readRDS("./inter_rds/02.k27ac.list.rds")

# check tissue order is the same! 
names(peaks.list)
names(k27ac.list)

# create function that returns subset based on hits
overlapped = function(query, queryhits.results){
  return(query[queryhits.results])
}

# function that returns peaks with NO overlapping hits
not.overlapping = function(query,queryhits.results){
  return(query[-queryhits.results])
}


# perform find overlaps, store results in list
overlap_results = purrr::map2(peaks.list, k27ac.list,~ queryHits(findOverlaps(.x,.y, type = "within")))

# subset based on queryHits and create new column
yes.k27ac = purrr::map2(peaks.list, overlap_results, ~ overlapped(.x,.y)) 
  
# add description column
yes.k27ac = purrr::map(yes.k27ac, ~ .x %>% add_mcol(.,"H3K27ac","+"))
  
# get non k27ac overlapping regions
no.k27ac = purrr::map2(peaks.list, overlap_results, ~ not.overlapping(.x, .y)) 
 
# add description column 
no.k27ac = purrr::map(no.k27ac, ~ .x %>% add_mcol(.,"H3K27ac","-"))
  
# combine granges list by tissue w/ annotated +/- K27ac
final = purrr::map2(yes.k27ac, no.k27ac, ~ append(.x,.y))
  
# create one granges
uniquePeaks.GR = c(final[[1]],final[[2]],final[[3]],final[[4]],final[[5]])

#-- annotate --# 
uniquePeaks.GR = annotatePeak(uniquePeaks.GR, 
                              TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                              annoDb = "org.Hs.eg.db")

#-- perform GO term analysis and plot --#
# get bg genes
hg19.genes = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

bg.genes = AnnotationDbi::select(org.Hs.eg.db, keytype = "ENTREZID", 
                                 columns = c("ENTREZID","SYMBOL"), 
                                 keys = unique(hg19.genes$gene_id))

# go analysis, data needs to be a dataframe
data = as.data.frame(uniquePeaks.GR@anno) %>%
  mutate(Tissue = factor(Tissue, levels = c("uniqueBM", "uniqueUC", "uniqueWAT", 
                                            "uniqueCH", "uniqueFB"))) %>%
  mutate(H3K27ac = factor(H3K27ac, levels = c("+","-")))

compareGO.tissue = compareCluster(SYMBOL~Tissue+H3K27ac,
                                  data = data,
                                  fun = "enrichGO", universe = bg.genes$SYMBOL,
                                  keyType = "SYMBOL", OrgDb = org.Hs.eg.db, 
                                  ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,
                                  readable = FALSE)

# plot
dotplot(compareGO.tissue, x = ~ H3K27ac, showCategory = 5) + 
  facet_grid(~ Tissue) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank()) +
  scale_radius(range = c(0.5,6)) +
  labs(title = "All Unique ATAC-seq Peaks", x = "H3K27ac")



