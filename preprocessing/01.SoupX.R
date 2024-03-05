setwd("/projects/b1042/LinLab/realign-mismatch")

samples<-paste0("EL-", 1:16)

### Run SoupX
library(tidyverse)
library(Seurat)
library(Signac)
library(SoupX)

for(i in samples){
  thissample=i
  
  ### Read filtered
  toc<-Read10X(paste0("/projects/b1042/LinLab/realign-mismatch/output/count/", thissample, "/outs/filtered_feature_bc_matrix/"))
  ### Read raw
  tod<-Read10X(paste0("/projects/b1042/LinLab/realign-mismatch/output/count/", thissample, "/outs/raw_feature_bc_matrix/"))
  #tod$`Gene Expression`
  #tod.split
  
  tod.seu<-CreateSeuratObject(tod$`Gene Expression`, min.cells = 0, min.features = 0)
  toc.seu<-CreateSeuratObject(toc$`Gene Expression`, min.cells = 0, min.features = 0)
  
  sc = SoupChannel(GetAssayData(tod.seu, assay = "RNA", slot = "counts"), GetAssayData(toc.seu, assay = "RNA", slot = "counts"))
  
  toc.seu %>% 
    SCTransform(vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = 0.7, verbose = FALSE) -> toc.seu
  
  DimPlot(toc.seu, reduction = "umap", label = T)+NoLegend()
  VlnPlot(toc.seu, features = c("nCount_RNA", "nFeature_RNA"))
  
  md <- toc.seu@meta.data
  sc = setClusters(sc, setNames(md$seurat_clusters, rownames(md)))
  sc = setDR(sc, toc.seu@reductions$umap@cell.embeddings)
  
  ### Calculate contamination fraction and adjust counts (roundToInt=T for SCTransform)
  #outpic=paste0("/projects/b1042/LinLab/realign-mismatch/output/SoupX/", thissample, "/SoupX-rho.png")
  #png(outpic, width = 800, height = 600)
  sc = autoEstCont(sc)
  #dev.off()
  
  out = adjustCounts(sc, roundToInt = T)
  
  cntSoggy = rowSums(sc$toc > 0)
  cntStrained = rowSums(out > 0)
  mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
  mostZeroed
  
  head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20)
  tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 20)
  
  #head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
  #plotMarkerDistribution(sc)
  
  ### Write decontaminated counts
  outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/SoupX/", thissample)
  dir.create(outname)
  DropletUtils::write10xCounts(path = outname, x = out, overwrite = T)
  
}

