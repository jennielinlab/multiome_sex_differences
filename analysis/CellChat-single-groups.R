library(tidyverse)
library(Seurat)
library(Signac)

library(CellChat)
library(future)

CellChatDB<-CellChatDB.mouse
showDatabaseCategory(CellChatDB)


library(NMF)
library(ggalluvial)


### Loop through groups
setwd("/projects/b1042/LinLab/realign-mismatch/scripts")

aggr.filt<-readRDS("../output/20000nCount/aggregated-20000nCountRNA-noMT-genes-noEndoPTcluster.rds")

groupsofinterest<-c("M.Young.Control", "M.Young.NTS", "F.Young.Control", "F.Young.NTS", "M.Aged.Control", "M.Aged.NTS", "F.Aged.Control", "F.Aged.NTS")


for(i in groupsofinterest){
  aggr.sub<-subset(aggr.filt, subset=sex.age.treat==i)
  md.sub<-aggr.sub@meta.data
  
  DefaultAssay(aggr.sub)<-"RNA"
  aggr.sub<-NormalizeData(aggr.sub)
  
  ### Get both RNA and SCT normalized data
  rna.norm<-aggr.sub@assays$RNA@data
  sct.norm<-aggr.sub@assays$SCT@data
  
  cc.obj<-createCellChat(object = rna.norm, meta = md.sub, group.by = "cellType2")
  
  cc.obj<-addMeta(cc.obj, md.sub)
  cc.obj<-setIdent(cc.obj, ident.use = "cellType2")
  groupSize<-as.numeric(table(cc.obj@idents))
  
  cc.obj@DB<-CellChatDB
  
  cc.obj<-subsetData(cc.obj)
  
  future::plan("multisession", workers = 12) # do parallel
  
  cc.obj <- identifyOverExpressedGenes(cc.obj)
  cc.obj <- identifyOverExpressedInteractions(cc.obj)
  
  # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  cc.obj <- projectData(cc.obj, PPI.mouse)
  
  ###
  cc.obj <- computeCommunProb(cc.obj)
  
  cc.obj.filt <- filterCommunication(cc.obj, min.cells = 10)
  
  ### Pathway level
  cc.obj.filt <- computeCommunProbPathway(cc.obj.filt)
  
  ### Calculate the aggregated cell-cell communication network
  cc.obj.filt <- aggregateNet(cc.obj.filt)
  
  groupSize <- as.numeric(table(cc.obj.filt@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cc.obj.filt@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cc.obj.filt@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  #
  
  # Compute the network centrality scores
  cc.obj.filt <- netAnalysis_computeCentrality(cc.obj.filt, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
  #netAnalysis_signalingRole_network(cc.obj.filt, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  
  outname=paste0("../output/20000nCount/CellChat/noMT/20000nCount-", i, ".rds")
  
  saveRDS(cc.obj.filt, file = outname)
  
  rm(aggr.sub)
  rm(cc.obj)
  rm(cc.obj.filt)
  gc()
  
  future::plan("sequential") # do parallel
  
}