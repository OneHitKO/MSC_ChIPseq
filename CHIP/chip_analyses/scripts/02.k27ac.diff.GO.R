#!/usr/bin/env Rscript

# This script focuses on the differential analysis results of H3K27ac 
# in BM vs FB/UC/WAT, CH vs FB/UC/WAT, and BM vs CH.
# Analyses performed: GO term analysis, signature identification
# Visualizations: volcano + violin plot, GO dot plot

library(tidyverse)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

#-- Data Import and Formatting--#
# import differential analysis results
k27ac.reports = readRDS("./rds/diffbind.k27ac.reports.rds")

# create list of reports focusing on BM contrasts or CH contrasts >> then data munge
k27ac_BM = k27ac.reports[str_which(names(k27ac.reports), pattern = "BM")]
k27ac_CH = k27ac.reports[-str_which(names(k27ac.reports), pattern = "BM")]

# only need fold change, genes, annotation, symbol, and new columns w/ Tissue
k27ac_BM = purrr::map(k27ac_BM, ~ .x %>% 
                        dplyr::mutate(Tissue = str_extract(names(.x)[[8]], pattern = "[:upper:]{2,3}")) %>%
                        dplyr::select("Fold","FDR","annotation","SYMBOL", "GENENAME","Tissue"))
  
k27ac_CH = purrr::map(k27ac_CH, ~ .x %>% 
                        dplyr::mutate(Tissue = str_extract(names(.x)[[8]], pattern = "[:upper:]{2,3}")) %>%
                        dplyr::select("Fold","FDR","annotation","SYMBOL", "GENENAME", "Tissue"))

# collapse into single data frame
k27ac_BM = do.call("rbind", k27ac_BM)
k27ac_CH = do.call("rbind", k27ac_CH)

# add CH vs BM by taking -1 * FC in k27ac_BM
CH = filter(k27ac_BM, Tissue == "CH") %>%
  dplyr::mutate(Fold = -1 * Fold, Tissue = "BM")

k27ac_CH = rbind(k27ac_CH, CH)

# create new simplified annotation column
simplify_order = c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Downstream", "Distal Intergenic")

k27ac_BM = k27ac_BM %>%
  dplyr::mutate(simplify = str_split_fixed(annotation, pattern = " \\(", n=2)[,1]) %>%
  dplyr::mutate(simplify = factor(simplify, levels = simplify_order))

k27ac_CH = k27ac_CH %>%
  dplyr::mutate(simplify = str_split_fixed(annotation, pattern = " \\(", n=2)[,1]) %>%
  dplyr::mutate(simplify = factor(simplify, levels = simplify_order))

#-- Visualization of Individual Peaks Data --#
# highlight these TFs
TFs = c("SOX9", "SOX11", "DLX5", "DLX6")

# violin BM with TFs highlighted
ggplot(k27ac_BM, aes(simplify, Fold))+
  geom_violin() +
  facet_grid(rows = vars(Tissue)) + 
  theme_bw() + 
  geom_dotplot(data = subset(k27ac_BM, SYMBOL %in% TFs),
               aes(fill = SYMBOL),
               binaxis = "y", stackdir = "center", dotsize = 1.5) +
  labs(title = "H3K27ac Differential Analysis: BM vs others")

# violin CH with TFs highlighted
ggplot(k27ac_CH, aes(simplify, Fold))+
  geom_violin() +
  facet_grid(rows = vars(Tissue)) + 
  theme_bw() + 
  geom_dotplot(data = subset(k27ac_CH, SYMBOL %in% TFs),
               aes(fill = SYMBOL),
               binaxis = "y", stackdir = "center", dotsize = 1.5) +
  labs(title = "H3K27ac Differential Analysis: CH vs others")


#-- Get Average H3K27ac FC For Every Genes --#
aveFold_k27ac_BM = k27ac_BM %>%
  dplyr::group_by(SYMBOL, GENENAME, Tissue) %>%
  dplyr::summarise(k27ac.av.Fold = mean(Fold))

aveFold_k27ac_CH = k27ac_CH %>%
  dplyr::group_by(SYMBOL, GENENAME, Tissue) %>%
  dplyr::summarise(k27ac.av.Fold = mean(Fold))

# Get genes that are always up in BM vs all others
up.inBM = aveFold_k27ac_BM %>%
  pivot_wider(names_from = Tissue, values_from = k27ac.av.Fold) %>%
  dplyr::filter(FB >=2, UC >=2, WAT >=2) %>%
  dplyr::mutate(Ave = (FB + UC + WAT)/3) %>%
  dplyr::arrange(desc(Ave)) %>%
  dplyr::select(SYMBOL, GENENAME, UC, WAT, FB, CH)

# Get genes that are always up in CH vs all others
up.inCH = aveFold_k27ac_CH %>%
  pivot_wider(names_from = Tissue, values_from = k27ac.av.Fold) %>%
  dplyr::filter(FB >=2, UC >=2, WAT >=2) %>%
  dplyr::mutate(Ave = (FB + UC + WAT)/3) %>%
  dplyr::arrange(desc(Ave)) %>%
  dplyr::select(SYMBOL, GENENAME, UC, WAT, FB, BM)

# Get common set of genes
up.inBM %>%
  dplyr::filter(SYMBOL %in% up.inCH$SYMBOL) %>%
  write_excel_csv("./reports/BM.CH.common.genes.csv")

# Get BM sig
up.inBM %>%
  dplyr::filter(!(SYMBOL %in% common.sig$SYMBOL)) %>%
  write_excel_csv("./reports/BM.only.genes.csv")

# Get CH sig
up.inCH %>%
  dplyr::filter(!(SYMBOL %in% common.sig$SYMBOL)) %>%
  write_excel_csv("./reports/CH.only.genes.csv")

#-- Perform GO Term Overrepresentation Analysis using average FC --#
# create new column describing if av FC is up (> 2) or down (< - 2)
aveFold_k27ac_BM = dplyr::mutate(aveFold_k27ac_BM,
                                 Direction = case_when(k27ac.av.Fold > 2 ~ "up in BM",
                                                       k27ac.av.Fold < -2 ~ "down in BM")) %>%
  dplyr::filter(!is.na(Direction))

aveFold_k27ac_CH = dplyr::mutate(aveFold_k27ac_CH, 
                                 Direction = case_when(k27ac.av.Fold > 2 ~ "up in CH",
                                                       k27ac.av.Fold < -2 ~ "down in CH")) %>%
  dplyr::filter(!is.na(Direction))

# change levels of factors
aveFold_k27ac_BM$Tissue = factor(aveFold_k27ac_BM$Tissue, levels = c("UC","WAT","FB","CH"))
aveFold_k27ac_BM$Direction = factor(aveFold_k27ac_BM$Direction, levels = c("up in BM", "down in BM"))

aveFold_k27ac_CH$Tissue = factor(aveFold_k27ac_CH$Tissue, levels = c("UC","WAT","FB"))
aveFold_k27ac_CH$Direction = factor(aveFold_k27ac_CH$Direction, levels = c("up in CH", "down in CH"))

# get bg genes for go analysis
hg19.genes = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
bg.genes = AnnotationDbi::select(org.Hs.eg.db, keytype = "ENTREZID", 
                                 columns = c("ENTREZID","SYMBOL"), 
                                 keys = unique(hg19.genes$gene_id))

# perform GO analysis for BM diff analysis, get csv
k27ac_compareGO_BM = compareCluster(SYMBOL~Tissue+Direction,
                                    data = aveFold_k27ac_BM, fun = "enrichGO", universe = bg.genes$SYMBOL,
                                    keyType = "SYMBOL", OrgDb = org.Hs.eg.db, 
                                    ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,
                                    readable = FALSE)

k27ac_GOresults_BM = k27ac_compareGO_BM@compareClusterResult

# perform GO analysis for CH diff analysis, get csv
k27ac_compareGO_CH = compareCluster(SYMBOL~Tissue+Direction,
                                    data = aveFold_k27ac_CH, fun = "enrichGO", universe = bg.genes$SYMBOL,
                                    keyType = "SYMBOL", OrgDb = org.Hs.eg.db, 
                                    ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05,
                                    readable = FALSE)

k27ac_GOresults_CH = k27ac_compareGO_CH@compareClusterResult


# plot top 8 results
dotplot(k27ac_compareGO_BM, x = ~Direction, showCategory = 8) + 
  facet_grid(~ Tissue) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank()) +
  scale_radius(range = c(0.5,6)) +
  labs(title = "GO Term Enrichment Analysis on Differential H3K27ac Peaks")

dotplot(k27ac_compareGO_CH, x = ~Direction, showCategory = 8) + 
  facet_grid(~ Tissue) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank()) +
  scale_radius(range = c(0.5,6)) +
  labs(title = "GO Term Enrichment Analysis on Differential H3K27ac Peaks")

# plot bone and cartilage related terms 
bone_cartilage_GOterms = unique(str_subset(k27ac_GOresults_BM$Description, 
                                    pattern = "osteo|bone|chondro|cartilage"))

# subset GO term enrichment analysis w/ just bone_cartilage terms
dplyr::filter(k27ac_GOresults_BM, Description %in% bone_cartilage_GOterms) %>%
  ggplot(., aes(Direction,Description))+
  geom_count(aes(color = qvalue, size = Count)) + 
  facet_grid(cols = vars(Tissue)) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        strip.background = element_rect(colour = "black")) +
  scale_color_viridis_c(direction = -1, option = "plasma", end = 0.85)
  
dplyr::filter(k27ac_GOresults_CH, Description %in% bone_cartilage_GOterms) %>%
  ggplot(., aes(Direction,Description))+
  geom_count(aes(color = qvalue, size = Count)) + 
  facet_grid(cols = vars(Tissue)) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        strip.background = element_rect(colour = "black")) +
  scale_color_viridis_c(direction = -1, option = "plasma", end = 0.85)
