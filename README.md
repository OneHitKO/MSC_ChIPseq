# Epigenetic profiling of different MSC cells
Can chromatin remodeling explain disaparities in chondro/osteogenic potential among different tissue-derived mesenchymal stromal cells?

Epigenetic profiles of undifferentiated  MSCs isolated from bone marrow, adipose, and umbilical
cord were compared with chondrocytes and fibroblasts. 

### Quick notes on software used for analysis
- Enrichment for histone post-translational modifications (H3K27ac, H3K4me1, H3K4me3) were called using `macs2` software on `.bam` files.

- Reproducible peaks for both ATAC-seq and ChIP-seq were called using `genrich`
software. 

- Differential analysis on ChIP-seq data was performed using `diffbind`. 

- Unique ATAC-seq peaks with or without overlapping reproducible H3K27ac peaks were identified using `bedtools`.

- Motif enrichment analysis was performed with `homer`.  
