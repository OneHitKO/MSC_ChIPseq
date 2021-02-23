#!/usr/bin/env Rscript

###--- Script will perform GO term analysis on motif family-specific regions ---###

library(tidyverse)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
library(clusterProfiler)

#-- import unique ATAC peaks --# 
dir = "../../peaks/uniquePeaks/"

bed.files = list.files(path=dir, pattern=".narrowPeak$", full.names=T)

# create list of granges of unique peaks
peaks.list = purrr::map(bed.files, ~rtracklayer::import(., format="narrowPeak"))

names(peaks.list) = str_extract(bed.files, pattern = "unique[:upper:]{2,3}")

#-- import motif files --# 
dir = "../../results_homer/"

# create nested list containing info on peaks enriched in motif families
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
# also just realized there was an easier output file I could've used in homer. (-m option)
# next time, redo the analysis w/ this homer command

GR.motifs = vector("list", length = 5)
names(GR.motifs) = names(motif.list)

for (i in seq_along(GR.motifs)){
  
  # initialize the nested list
  GR.motifs[[i]] = vector("list", length(motif.list[[i]]))
  names(GR.motifs[[i]]) = names(motif.list[[i]])

  # subset granges, probably a better way to do this
  # note: try map2, see below when submitting on k27ac 
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

names(k27ac.list)=str_extract(k27ac.files, pattern ="[:upper:]{2,3}_K27ac")

# create functions to use w/ map or lapply
# add a new column to metadata of granges
add_mcol = function(gr_obj,col_name,value){
  mcols(gr_obj)[[toString(col_name)]] = value
  return(gr_obj)
}

# return based on hits
overlapped = function(query, queryhits.results){
  return(query[queryhits.results])
}

# return NOT hits
not.overlapping = function(query,queryhits.results){
  return(query[-queryhits.results])
}

# initialize lists
overlap_results = vector("list", length = 5)
yes.k27ac = vector("list", length = 5)
no.k27ac = vector("list", length = 5)
final = vector("list", length = 5)
names(final) = names(GR.motifs)

for (i in seq_along(GR.motifs)){

  # find overlaps between motif regions and k27ac peaks
  overlap_results[[i]] = purrr::map2(GR.motifs[[i]], k27ac.list[i],
                                ~ queryHits(findOverlaps(.x,.y, type = "within")))
  
  # subset based on queryHits and create new column
  yes.k27ac[[i]] = purrr::map2(GR.motifs[[i]], overlap_results[[i]], ~ overlapped(.x, .y)) 
  
  yes.k27ac[[i]] = purrr::map(yes.k27ac[[i]], ~ .x %>% 
                                add_mcol(.,"H3K27ac","+"))
  
  # get non k27ac overlapping regions
  no.k27ac[[i]] = purrr::map2(GR.motifs[[i]], overlap_results[[i]], ~ not.overlapping(.x, .y)) 
  
  no.k27ac[[i]] = purrr::map(no.k27ac[[i]], ~ .x %>% 
                               add_mcol(.,"H3K27ac","-"))
  
  # combine granges list by tissue w/ annotated +/- K27ac
  final[[i]] = purrr::map2(yes.k27ac[[i]], no.k27ac[[i]], ~ append(.x,.y))
  
  # add additinal metadata columns based on motif, family, tissue
  final[[i]] = purrr::map(final[[i]], ~ .x %>%
                            add_mcol(., "Tissue",names(final)[[i]]))
  
}

#-- annotate --# 

# function for nested map
annotate = function(x){
  purrr::map(x, ~ annotatePeak(peak = .x, 
                               TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                               annoDb = "org.Hs.eg.db"))
} 

final = purrr::map(final, annotate)

#-- create individual lists to prepare for GO term analysis --#
bzip = c(final$uniqueBM$`Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer(0.990)`@anno,
         final$uniqueCH$`AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer(0.962)`@anno,
         final$uniqueFB$`BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer(0.995)`@anno,
         final$uniqueUC$`Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer(0.992)`@anno,
         final$uniqueWAT$`BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer(0.996)`@anno) %>%
  as.data.frame() %>%
  mutate(H3K27ac = factor(H3K27ac, levels = c("+", "-")))

runt = c(final$uniqueBM$`RUNX2(Runt)/PCa-RUNX2-ChIP-Seq(GSE33889)/Homer(0.942)`@anno,
         final$uniqueWAT$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer(0.967)`@anno) %>%
  as.data.frame() %>%
  mutate(H3K27ac = factor(H3K27ac, levels = c("+", "-")))

tead = c(final$uniqueBM$`TEAD(TEA)/Fibroblast-PU.1-ChIP-Seq(Unpublished)/Homer(0.971)`@anno,
         final$uniqueUC$`TEAD1(TEAD)/HepG2-TEAD1-ChIP-Seq(Encode)/Homer(0.914)`@anno,
         final$uniqueWAT$`TEAD2/MA1121.1/Jaspar(0.950)`@anno) %>%
  as.data.frame() %>%
  mutate(H3K27ac = factor(H3K27ac, levels = c("+", "-")))

#-- perform GO term analysis and plot --#
# get bg genes
hg19.genes = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

bg.genes = AnnotationDbi::select(org.Hs.eg.db, keytype = "ENTREZID", 
                                 columns = c("ENTREZID","SYMBOL"), 
                                 keys = unique(hg19.genes$gene_id))

# bzip motif 
compareGO.bzip = compareCluster(SYMBOL~Tissue+H3K27ac,
                           data = bzip, fun = "enrichGO", universe = bg.genes$SYMBOL,
                           keyType = "SYMBOL", OrgDb = org.Hs.eg.db, 
                           ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,
                           readable = FALSE)

dotplot(compareGO.bzip, x = ~ H3K27ac, showCategory = 5) + 
  facet_grid(~ Tissue) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank()) +
  scale_radius(range = c(0.5,6)) +
  labs(x = "H3K27ac", title = "bZIP motif")

# runt motif
compareGO.runt = compareCluster(SYMBOL~Tissue+H3K27ac,
                                data = runt, fun = "enrichGO", universe = bg.genes$SYMBOL,
                                keyType = "SYMBOL", OrgDb = org.Hs.eg.db, 
                                ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,
                                readable = FALSE)

dotplot(compareGO.runt, x = ~ H3K27ac, showCategory = 8) + 
  facet_grid(~ Tissue) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank()) +
  scale_radius(range = c(0.5,6)) +
  labs(x = "H3K27ac", title = "Runt motif")

# tead motif
compareGO.tead = compareCluster(SYMBOL~Tissue+H3K27ac,
                                data = tead, fun = "enrichGO", universe = bg.genes$SYMBOL,
                                keyType = "SYMBOL", OrgDb = org.Hs.eg.db, 
                                ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,
                                readable = FALSE)

dotplot(compareGO.tead, x = ~ H3K27ac, showCategory = 8) + 
  facet_grid(~ Tissue) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank()) +
  scale_radius(range = c(0.5,6)) +
  labs(x = "H3K27ac", title = "TEA motif")



