# Get DEGs between groups + celltypes

library(tidyverse)
library(Seurat)
library(Signac)

### Read aggregated harmony-corrected seurat object
aggr.harm<-readRDS("/projects/b1042/LinLab/realign-mismatch/output/aggr/aggregated-seurat-SoupX-20000nCountRNA-dedoub-harmony.rds")

Idents(aggr.harm)<-"cellType"

DefaultAssay(aggr.harm)<-"RNA"
aggr.harm<-NormalizeData(aggr.harm)

for(i in unique(Idents(aggr.harm))){
  celltype=as.character(i)
  
  ### Male vs Female
  ### Young male control vs. young female control
  if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Young.Control"]==i)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Young.Control"]==i)>=3){
    de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Young.Control", ident.2 = "M.Young.Control", subset.ident = celltype, test.use = "LR", assay = "RNA")
    de.res %>% rownames_to_column("Gene") -> de.res
    
    outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/genes/MvF/Young/Control/", celltype, "-young.F.control-vs-young.M.control.txt")
    write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
  }
  
  ### Young male NTS vs. young female NTS
  if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Young.NTS"]==i)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Young.NTS"]==i)>=3){
    de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Young.NTS", ident.2 = "M.Young.NTS", subset.ident = celltype, test.use = "LR", assay = "RNA")
    de.res %>% rownames_to_column("Gene") -> de.res
    
    outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/genes/MvF/Young/NTS/", celltype, "-young.F.NTS-vs-young.M.NTS.txt")
    write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
  }
  
  ### Aged male control vs. aged female control
  if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Aged.Control"]==i)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Aged.Control"]==i)>=3){
    de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Aged.Control", ident.2 = "M.Aged.Control", subset.ident = celltype, test.use = "LR", assay = "RNA")
    de.res %>% rownames_to_column("Gene") -> de.res
    
    outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/genes/MvF/Aged/Control/", celltype, "-aged.F.control-vs-aged.M.control.txt")
    write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
  }
  
  ### Aged male NTS vs. aged female NTS
  if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Aged.NTS"]==i)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Aged.NTS"]==i)>=3){
    de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Aged.NTS", ident.2 = "M.Aged.NTS", subset.ident = celltype, test.use = "LR", assay = "RNA")
    de.res %>% rownames_to_column("Gene") -> de.res
    
    outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/genes/MvF/Aged/NTS/", celltype, "-aged.F.NTS-vs-aged.M.NTS.txt")
    write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
  }
  
  
  ### Young vs. aged
  ### Young female control vs. Aged female control
  if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Young.Control"]==i)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Aged.Control"]==i)>=3){
    de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Young.Control", ident.2 = "F.Aged.Control", subset.ident = celltype, test.use = "LR", assay = "RNA")
    de.res %>% rownames_to_column("Gene") -> de.res
    
    outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/genes/YvO/F/Control/", celltype, "-young.F.control-vs-aged.F.control.txt")
    write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
  }
  ### Young female NTS vs. Aged female NTS
  if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Young.NTS"]==i)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Aged.NTS"]==i)>=3){
    de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Young.NTS", ident.2 = "F.Aged.NTS", subset.ident = celltype, test.use = "LR", assay = "RNA")
    de.res %>% rownames_to_column("Gene") -> de.res
    
    outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/genes/YvO/F/NTS/", celltype, "-young.F.NTS-vs-aged.F.NTS.txt")
    write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
  }
  
  ### Young male control vs. Aged male control
  if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Young.Control"]==i)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Aged.Control"]==i)>=3){
    de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "M.Young.Control", ident.2 = "M.Aged.Control", subset.ident = celltype, test.use = "LR", assay = "RNA")
    de.res %>% rownames_to_column("Gene") -> de.res
    
    outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/genes/YvO/M/Control/", celltype, "-young.M.control-vs-aged.M.control.txt")
    write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
  }
  ### Young male NTS vs. Aged male NTS
  if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Young.NTS"]==i)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Aged.NTS"]==i)>=3){
    de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "M.Young.NTS", ident.2 = "M.Aged.NTS", subset.ident = celltype, test.use = "LR", assay = "RNA")
    de.res %>% rownames_to_column("Gene") -> de.res
    
    outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/genes/YvO/M/NTS/", celltype, "-young.M.NTS-vs-aged.M.NTS.txt")
    write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
  }
  
  ### Control vs. NTS
  ### Young female control vs. young female NTS
  if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Young.Control"]==i)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Young.NTS"]==i)>=3){
    de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Young.Control", ident.2 = "F.Young.NTS", subset.ident = celltype, test.use = "LR", assay = "RNA")
    de.res %>% rownames_to_column("Gene") -> de.res
    
    outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/genes/CvN/F/Young/", celltype, "-young.F.control-vs-young.F.NTS.txt")
    write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
  }
  ### Aged female control vs. Aged female NTS
  if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Aged.Control"]==i)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Aged.NTS"]==i)>=3){
    de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Aged.Control", ident.2 = "F.Aged.NTS", subset.ident = celltype, test.use = "LR", assay = "RNA")
    de.res %>% rownames_to_column("Gene") -> de.res
    
    outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/genes/CvN/F/Aged/", celltype, "-aged.F.control-vs-aged.F.NTS.txt")
    write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
  }
  
  ### Young male control vs. young male NTS
  if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Young.Control"]==i)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Young.NTS"]==i)>=3){
    de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "M.Young.Control", ident.2 = "M.Young.NTS", subset.ident = celltype, test.use = "LR", assay = "RNA")
    de.res %>% rownames_to_column("Gene") -> de.res
    
    outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/genes/CvN/M/Young/", celltype, "-young.M.control-vs-young.M.NTS.txt")
    write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
  }
  ### Aged male control vs. Aged male NTS
  if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Aged.Control"]==i)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Aged.NTS"]==i)>=3){
    de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "M.Aged.Control", ident.2 = "M.Aged.NTS", subset.ident = celltype, test.use = "LR", assay = "RNA")
    de.res %>% rownames_to_column("Gene") -> de.res
    
    outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/genes/CvN/M/Aged/", celltype, "-aged.M.control-vs-aged.M.NTS.txt")
    write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
  }
  
}

### Get celltype markers
allmarkers<-FindAllMarkers(object = aggr.harm, test.use = "LR")

write.table(allmarkers, "/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/genes/celltype/AllMarkers-RNA.txt", sep="\t", col.names=T, row.names=F, quote=F)

