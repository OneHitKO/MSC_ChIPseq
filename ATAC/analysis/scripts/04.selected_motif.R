#!/bin/bash/ Rscript

# This script takes the results of a motif enrichment analysis which identified
# instances of motif of interest in all ATAC-seq peaks stratified by
# overlap with H3K27ac.
# 1. Annotate all reproducible ATAC-seq peaks from each cell-type with Homer 
#   - Motifs of interest: ATF3, DLX1, DLX5, FOSL2, JUN, RUNX1, RUNX2, SOX9
# 2. In R: for each motif, create list with output of each tissue. Filter out
#    peaks with NO instance of motif. 
# 3. Find regions with and without overlap to H3K27ac, perform GO term
#    analysis. 

library(tidyverse)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
library(clusterProfiler)
library(shiny)

##-- Data import, filtering, and reformatting --##
# get directory of the homer output text files 
dir = "~/work/msc/ATAC/results_homer/skeletal_motif_regions/"
files = list.files(path = dir, full.names = T)

# create list of all files for filtering
list.all = map(files, read_tsv)

# rename files 
names(list.all) = basename(files)

# remove peaks that do not contain motif of interest,
# which the info is in the last column (22)
# also add new column with tissue
list.all = map(list.all, ~.x %>%
               dplyr::rename(Motif_Instance = 22, Peak_ID = 1) %>%
               dplyr::filter(!is.na(Motif_Instance)))

# convert to granges
list.all = map(list.all, ~ makeGRangesFromDataFrame(.x, keep.extra.columns = T))


##-- Overlapping with H3K27ac peaks --##
# import H3K27ac peaks data
k27ac.list = readRDS("./inter_rds/02.k27ac.list.rds")

# try overlapping without nested lists and map2 :P 
# this code is so much shorter and more simplified! should update others.
for (i in seq_along(list.all)){
  # define motif GR as query
  query.motif = list.all[[i]]
  
  # get tissue name
  tissue = word(names(list.all)[i], start = 1, sep = "\\.")
  
  # get TF name
  tf = word(names(list.all)[i], start = 2, sep = "\\.")
  
  # define appropriate k27ac GR as subject
  subject.k27ac = k27ac.list[[str_which(names(k27ac.list), pattern = tissue)]]
  
  # perform overlaps
  hits = findOverlaps(query.motif, subject.k27ac, type = "within")
  
  # initiate new H3K27ac column
  mcols(list.all[[i]])$H3K27ac = NA
  
  # subset GRanges based on overlapping hits, add new mcol for H3K27ac
  mcols(list.all[[i]][queryHits(hits)])$H3K27ac = "+"
  mcols(list.all[[i]][-queryHits(hits)])$H3K27ac = "-"
  
  # add tissue column
  mcols(list.all[[i]])$Tissue = tissue
  
  # add motif column
  mcols(list.all[[i]])$TF = tf
}

# convert granges to df
list.all = map(list.all, as.data.frame)

# create one large data frame
final.df = do.call(rbind, list.all)

# reset levels
final.df = final.df %>%
  dplyr::mutate(H3K27ac = factor(H3K27ac, levels = c("+","-")),
                Tissue = factor(Tissue, levels = c("BM","CH","UC","WAT","FB")))


##-- GO term analysis --##
# initiate list for go term enrichment results, ordred by TF
tfs = unique(final.df$TF)
go.results = vector("list", length(tfs))
names(go.results) = tfs

# get list of bg genes for analysis
hg19 = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
bg.genes = AnnotationDbi::select(org.Hs.eg.db, keytype = "ENTREZID", 
                                 columns = c("ENTREZID","SYMBOL"), 
                                 keys = unique(hg19$gene_id))

# perform go term analysis
for (i in seq_along(go.results)){
  # first subset df with motif results/annotation for tf of interest
  go.results[[i]] = final.df %>%
    dplyr::filter(TF == names(go.results)[i]) %>%
    compareCluster(Gene.Name~Tissue+H3K27ac,
                   data = ., fun = "enrichGO", universe = bg.genes$SYMBOL,
                   keyType = "SYMBOL", OrgDb = org.Hs.eg.db, 
                   ont = "BP", pAdjustMethod = "BH", readable = FALSE)
}

##-- Plots --##
# create shiny app to visualize all results
# define UI 
ui = fluidPage(
  sidebarLayout(
    sidebarPanel(
      
      # select TF of interest
      selectInput(inputId = "tfactor", label = "Select transcription factor of interest",
                  choices = names(go.results))
    ),
    
    mainPanel(plotOutput(outputId = "dotplot"))
  )
)

# server logic
server = function(input, output){
  
  output$dotplot = renderPlot({
    # plot dot plot of tf of interest of top 5 terms for each group
    dotplot(go.results[[str_which(names(go.results), pattern = input$tfactor)]],
            x = ~ H3K27ac,
            showCategory = 5) +
      facet_grid(~ Tissue) +
      theme(axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 10),
            panel.grid = element_blank()) +
      labs(x = "H3K27ac", title = input$tfactor)
  })
}

# run app
shinyApp(ui = ui, server = server)





