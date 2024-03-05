setwd("/projects/b1042/LinLab/realign-mismatch")

library(tidyverse)
library(Seurat)
library(Signac)
library(DoubletFinder)

### DoubletFinder loop
### Adapted from Muto et al. 2021 (https://github.com/p4rkerw/Muto_Wilson_NComm_2020/blob/master/snRNA_prep/seurat_rna_process.R)
doubletFunction<-function(seu.obj, libraryID){
  seu.sub <- subset(seu.obj, subset=sample==libraryID)
  seu.sub <- SCTransform(seu.sub, vst.flavor = "v2", verbose = FALSE)
  seu.sub <- RunPCA(seu.sub, verbose = FALSE) 
  seu.sub <- RunUMAP(seu.sub, dims = 1:20)
  
  ### Expected doublet rate = 0.8% per 1000 cells (10x)
  exp.doublet.rate <- (0.008/1000) * ncol(seu.sub)
  
  sweep.res.list <- paramSweep_v3(seu.sub, PCs = 1:20, sct = T)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # select the pK that corresponds to max bcmvn to optimize doublet detection
  pK <- bcmvn %>% 
    filter(BCmetric == max(BCmetric)) %>% 
    dplyr::select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  
  nExp_poi <- round(exp.doublet.rate*ncol(seu.sub))  ## Doublet formation rate based on 0.8% doublets per 1000 cells. 
  
  seu.sub <- doubletFinder_v3(seu.sub, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
  doub.md<-seu.sub@meta.data %>% rownames_to_column("CellID")
  
  doublet.out<-data.frame(sample=doub.md$CellID,
                          doublet.ann=doub.md[,ncol(doub.md)])
  
  return(doublet.out)
}



aggr.filt2<-readRDS("/projects/b1042/LinLab/realign-mismatch/output/aggr/aggregated-seurat-SoupX-relaxedCutoff-20000nCountRNA.rds")

### Already filtered
#aggr.filt<-subset(aggr.seu, 
#                  subset=percent.mt<1 & 
#                    #percent.ribo<5 & 
#                    nCount_RNA>350 & 
#                    nCount_RNA<20000 &
#                    nFeature_RNA > 200 &
#                    #nFeature_RNA < 4000 & 
#                    nCount_ATAC > 1000 &
#                    #nCount_ATAC < 100000 &
#                    nucleosome_signal < 2 &
#                    TSS.enrichment > 2)

ncol(aggr.filt2)
table(aggr.filt2$sample)
### max(table(aggr.filt$sample)) - max, 7025 cells
### min(table(aggr.filt$sample)) - min, 910 cells

rm(aggr.seu.replaced)
gc()

#samples<-c("EL-", 1:16)

sampleids<-unique(aggr.filt2$sample)


doublet.anns<-lapply(sampleids, FUN = function(x) {doubletFunction(seu.obj = aggr.filt2, libraryID = x)})

doublet.ann.df<-bind_rows(doublet.anns)

write.table(doublet.ann.df, "/projects/b1042/LinLab/realign-mismatch/output/doublets-5percent-relaxed-noRiboFilter-20000nCountRNA.txt", sep="\t", col.names = T, row.names = F, quote=F)
