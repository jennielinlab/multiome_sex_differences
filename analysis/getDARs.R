# Get DAR between groups + celltypes

library(tidyverse)
library(Seurat)
library(Signac)

args = commandArgs(trailingOnly=TRUE)


### Read aggregated harmony-corrected seurat object
aggr.harm<-readRDS("/projects/b1042/LinLab/realign-mismatch/output/aggr/aggregated-seurat-SoupX-20000nCountRNA-dedoub-harmony-MACS-noATAC-celltype.rds")

Idents(aggr.harm)<-"cellType"

DefaultAssay(aggr.harm)<-"peaks"

celltype=args[1]

### Male vs Female
### Young male control vs. young female control
if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Young.Control"]==celltype)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Young.Control"]==celltype)>=3){
  de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Young.Control", ident.2 = "M.Young.Control", subset.ident = celltype, test.use = "LR", assay = "peaks", min.pct = 0.05, latent.vars = 'nCount_peaks')
  de.res %>% rownames_to_column("Peak") -> de.res
  
  outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/peaks/MvF/Young/Control/", celltype, "-young.F.control-vs-young.M.control.txt")
  write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
}

### Young male NTS vs. young female NTS
if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Young.NTS"]==celltype)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Young.NTS"]==celltype)>=3){
  de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Young.NTS", ident.2 = "M.Young.NTS", subset.ident = celltype, test.use = "LR", assay = "peaks", min.pct = 0.05, latent.vars = 'nCount_peaks')
  de.res %>% rownames_to_column("Peak") -> de.res
  
  outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/peaks/MvF/Young/NTS/", celltype, "-young.F.NTS-vs-young.M.NTS.txt")
  write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
}

### Aged male control vs. aged female control
if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Aged.Control"]==celltype)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Aged.Control"]==celltype)>=3){
  de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Aged.Control", ident.2 = "M.Aged.Control", subset.ident = celltype, test.use = "LR", assay = "peaks", min.pct = 0.05, latent.vars = 'nCount_peaks')
  de.res %>% rownames_to_column("Peak") -> de.res
  
  outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/peaks/MvF/Aged/Control/", celltype, "-aged.F.control-vs-aged.M.control.txt")
  write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
}

### Aged male NTS vs. aged female NTS
if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Aged.NTS"]==celltype)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Aged.NTS"]==celltype)>=3){
  de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Aged.NTS", ident.2 = "M.Aged.NTS", subset.ident = celltype, test.use = "LR", assay = "peaks", min.pct = 0.05, latent.vars = 'nCount_peaks')
  de.res %>% rownames_to_column("Peak") -> de.res
  
  outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/peaks/MvF/Aged/NTS/", celltype, "-aged.F.NTS-vs-aged.M.NTS.txt")
  write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
}


### Young vs. aged
### Young female control vs. Aged female control
if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Young.Control"]==celltype)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Aged.Control"]==celltype)>=3){
  de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Young.Control", ident.2 = "F.Aged.Control", subset.ident = celltype, test.use = "LR", assay = "peaks", min.pct = 0.05, latent.vars = 'nCount_peaks')
  de.res %>% rownames_to_column("Peak") -> de.res
  
  outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/peaks/YvO/F/Control/", celltype, "-young.F.control-vs-aged.F.control.txt")
  write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
}
### Young female NTS vs. Aged female NTS
if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Young.NTS"]==celltype)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Aged.NTS"]==celltype)>=3){
  de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Young.NTS", ident.2 = "F.Aged.NTS", subset.ident = celltype, test.use = "LR", assay = "peaks", min.pct = 0.05, latent.vars = 'nCount_peaks')
  de.res %>% rownames_to_column("Peak") -> de.res
  
  outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/peaks/YvO/F/NTS/", celltype, "-young.F.NTS-vs-aged.F.NTS.txt")
  write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
}

### Young male control vs. Aged male control
if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Young.Control"]==celltype)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Aged.Control"]==celltype)>=3){
  de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "M.Young.Control", ident.2 = "M.Aged.Control", subset.ident = celltype, test.use = "LR", assay = "peaks", min.pct = 0.05, latent.vars = 'nCount_peaks')
  de.res %>% rownames_to_column("Peak") -> de.res
  
  outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/peaks/YvO/M/Control/", celltype, "-young.M.control-vs-aged.M.control.txt")
  write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
}
### Young male NTS vs. Aged male NTS
if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Young.NTS"]==celltype)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Aged.NTS"]==celltype)>=3){
  de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "M.Young.NTS", ident.2 = "M.Aged.NTS", subset.ident = celltype, test.use = "LR", assay = "peaks", min.pct = 0.05, latent.vars = 'nCount_peaks')
  de.res %>% rownames_to_column("Peak") -> de.res
  
  outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/peaks/YvO/M/NTS/", celltype, "-young.M.NTS-vs-aged.M.NTS.txt")
  write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
}

### Control vs. NTS
### Young female control vs. young female NTS
if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Young.Control"]==celltype)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Young.NTS"]==celltype)>=3){
  de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Young.Control", ident.2 = "F.Young.NTS", subset.ident = celltype, test.use = "LR", assay = "peaks", min.pct = 0.05, latent.vars = 'nCount_peaks')
  de.res %>% rownames_to_column("Peak") -> de.res
  
  outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/peaks/CvN/F/Young/", celltype, "-young.F.control-vs-young.F.NTS.txt")
  write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
}
### Aged female control vs. Aged female NTS
if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Aged.Control"]==celltype)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="F.Aged.NTS"]==celltype)>=3){
  de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "F.Aged.Control", ident.2 = "F.Aged.NTS", subset.ident = celltype, test.use = "LR", assay = "peaks", min.pct = 0.05, latent.vars = 'nCount_peaks')
  de.res %>% rownames_to_column("Peak") -> de.res
  
  outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/peaks/CvN/F/Aged/", celltype, "-aged.F.control-vs-aged.F.NTS.txt")
  write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
}

### Young male control vs. young male NTS
if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Young.Control"]==celltype)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Young.NTS"]==celltype)>=3){
  de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "M.Young.Control", ident.2 = "M.Young.NTS", subset.ident = celltype, test.use = "LR", assay = "peaks", min.pct = 0.05, latent.vars = 'nCount_peaks')
  de.res %>% rownames_to_column("Peak") -> de.res
  
  outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/peaks/CvN/M/Young/", celltype, "-young.M.control-vs-young.M.NTS.txt")
  write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
}
### Aged male control vs. Aged male NTS
if(sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Aged.Control"]==celltype)>=3 & sum(aggr.harm$cellType[aggr.harm$sex.age.treat=="M.Aged.NTS"]==celltype)>=3){
  de.res<-FindMarkers(aggr.harm, group.by = "sex.age.treat", ident.1 = "M.Aged.Control", ident.2 = "M.Aged.NTS", subset.ident = celltype, test.use = "LR", assay = "peaks", min.pct = 0.05, latent.vars = 'nCount_peaks')
  de.res %>% rownames_to_column("Peak") -> de.res
  
  outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/peaks/CvN/M/Aged/", celltype, "-aged.M.control-vs-aged.M.NTS.txt")
  write.table(de.res, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
}

### Get celltype markers

allmarkers<-FindMarkers(object = aggr.harm, ident.1 = celltype, test.use = "LR", assay = "peaks", min.pct = 0.05, latent.vars = 'nCount_peaks')
allmarkers %>% rownames_to_column("Peak") -> allmarkers

outname=paste0("/projects/b1042/LinLab/realign-mismatch/output/DE/20000nCountRNA/peaks/celltype/", celltype, "-peaks.txt")
write.table(allmarkers, file = outname, sep="\t", col.names=T, row.names=F, quote=F)

