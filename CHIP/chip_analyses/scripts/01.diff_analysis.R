#!/usr/bin/env Rscript

## This script will perform differential enrichment analysis of ChIP-seq peaks
## of bone marrow-derived MSCs (BM-MSC) vs other MSCs and cell types.

library("DiffBind")
#library("tidyverse")

#-- read in peaksets  
k27ac = dba(sampleSheet="diffbind.k27ac.csv")
#k4m3 = dba(sampleSheet="diffbind.k4m3.csv")
#k4m1 = dba(sampleSheet="diffbind.k4m1.csv")

#-- calculate binding matrix with consensus peakset
k27ac = dba.count(k27ac, bParallel=FALSE)
#k4m3 = dba.count(k4m3)
#k4m1 = dba.count(k4m1)

#-- data normalization by library size, without conservative background subtr.
k27ac = dba.normalize(k27ac, normalize=DBA_NORM_LIB)
#k4m3 = dba.normalize(k4m3, normalize=DBA_NORM_LIB)
#k4m1 = dba.normalize(k4m1, normalize=DBA_NORM_LIB)

#-- contrasts for k27ac
# contrast1: BM vs oth
BM.vs.FB = dba.contrast(k27ac, group1 = k27ac$masks$BM, 
                        group2 = k27ac$masks$FB, 
                        name1 = "BM", name2 = "FB")

BM.vs.UC = dba.contrast(k27ac, group1 = k27ac$masks$BM,
                        group2 = k27ac$masks$UC,
                        name1 = "BM", name2 = "UC")

BM.vs.WAT = dba.contrast(k27ac, group1 = k27ac$masks$BM,
                         group2 = k27ac$masks$WAT,
                         name1 = "BM", name2 = "WAT")

# contrast2: CH vs oth
CH.vs.FB = dba.contrast(k27ac, group1 = k27ac$masks$CH,
                        group2 = k27ac$masks$FB,
                        name1 = "CH", name2 = "FB")

CH.vs.UC = dba.contrast(k27ac, group1 = k27ac$masks$CH,
                        group2 = k27ac$masks$UC,
                        name1 = "CH", name2 = "UC")

CH.vs.WAT = dba.contrast(k27ac, group1 = k27ac$masks$CH,
                         group2 = k27ac$masks$WAT,
                         name1 = "CH", name2 = "WAT")

# contrast3: BM vs CH
BM.vs.CH = dba.contrast(k27ac, group1 = k27ac$masks$BM,
                         group2 = k27ac$masks$CH,
                         name1 = "BM", name2 ="CH")

# create list of contrasts
contrast.list = list(BM.vs.FB, BM.vs.UC, BM.vs.WAT, 
                     CH.vs.FB, CH.vs.UC, CH.vs.WAT, BM.vs.CH)

names(contrast.list) = c("BM.vs.FB", "BM.vs.UC", "BM.vs.WAT", 
                         "CH.vs.FB", "CH.vs.UC", "CH.vs.WAT", "BM.vs.CH")

#-- analyze k27ac
k27ac.contrast.list = purrr::map(contrast.list, ~ dba.analyze(.x))


