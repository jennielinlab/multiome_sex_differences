###

library(tidyverse)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(GenomicFeatures)
library(EnsDb.Mmusculus.v79)

#annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
#seqlevelsStyle(annotation) <- "UCSC"
annotation<-readRDS("/projects/b1042/LinLab/realign-mismatch/data/EnsDb.Mmusculus.v79-annotation.rds")


### Read in aggregated samples
aggr<-Read10X("../output/aggr/NTS_combined/outs/filtered_feature_bc_matrix/")


aggr.seu<-CreateSeuratObject(counts = aggr$`Gene Expression`, 
                             min.cells = 0, min.features = 0, assay = "RNA")
#aggr.seu<-CreateSeuratObject(counts = aggr$`Gene Expression`, 
#                              min.cells = 3, min.features = 200, assay = "RNA")

fragpath<-paste0("/projects/b1042/LinLab/realign-mismatch/output/aggr/NTS_combined/outs/atac_fragments.tsv.gz")

# create ATAC assay and add it to the object
aggr.seu[["ATAC"]] <- CreateChromatinAssay(
  counts = aggr$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

aggr.seu[["percent.mt"]] <- PercentageFeatureSet(aggr.seu, pattern = "^mt-")
aggr.seu[["percent.rps"]] <- PercentageFeatureSet(aggr.seu, pattern = "^Rps")
aggr.seu[["percent.rpl"]] <- PercentageFeatureSet(aggr.seu, pattern = "^Rpl")
aggr.seu$percent.ribo<-aggr.seu$percent.rps+aggr.seu$percent.rpl


aggr.seu@meta.data %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt))+
  geom_point()+
  scale_color_viridis_c()
ggsave("../plots/QC/All-samples-nCount_RNA-vs-nFeature_RNA.png", width = 8, height = 6)

aggr.seu@meta.data %>%
  ggplot(aes(x=nCount_ATAC, y=nFeature_ATAC, color=percent.mt))+
  geom_point()+
  scale_color_viridis_c()
ggsave("../plots/QC/All-samples-nCount_ATAC-vs-nFeature_ATAC.png", width = 8, height = 6)

VlnPlot(
  object = aggr.seu,
  features = c("nCount_RNA", "nCount_ATAC", "percent.mt", "percent.ribo"),
  ncol = 2,
  pt.size = 0.2
)
ggsave("../plots/QC/All-samples-violin-plots.png", width = 8, height = 6)

DefaultAssay(aggr.seu) <- "ATAC"

aggr.seu <- NucleosomeSignal(aggr.seu)
aggr.seu <- TSSEnrichment(aggr.seu)

VlnPlot(
  object = aggr.seu,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt", "percent.ribo"),
  ncol = 3,
  pt.size = 0
)
ggsave("../plots/QC/All-samples-violin-plots-v2.png", width = 8, height = 6)

aggr.seu$sample<-unlist(lapply(str_split(colnames(aggr.seu), pattern = "-"), "[[", 2))
table(aggr.seu$sample)

VlnPlot(object = aggr.seu, group.by = "sample", features = "percent.mt")

#saveRDS(aggr.seu, "/projects/b1042/LinLab/realign-mismatch/output/aggr/aggregated-seurat.rds")

aggr.seu<-readRDS("/projects/b1042/LinLab/realign-mismatch/output/aggr/aggregated-seurat.rds")

VlnPlot(
  object = aggr.seu,
  features = c("percent.mt"),
  ncol = 3,
  pt.size = 0, y.max = 5
)

### Run SoupX
# SoupX.R

### Read in all SoupX counts, merge, and replace counts in aggregated matrix
el1<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-1/")
el1.seu<-CreateSeuratObject(el1)
el1.seu<-RenameCells(el1.seu, new.names = paste0(unlist(lapply(str_split(colnames(el1.seu), pattern = "-"), "[[", 1)), -1))

el2<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-2/")
el2.seu<-CreateSeuratObject(el2)
el2.seu<-RenameCells(el2.seu, new.names = paste0(unlist(lapply(str_split(colnames(el2.seu), pattern = "-"), "[[", 1)), -2))

el3<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-3/")
el3.seu<-CreateSeuratObject(el3)
el3.seu<-RenameCells(el3.seu, new.names = paste0(unlist(lapply(str_split(colnames(el3.seu), pattern = "-"), "[[", 1)), -3))

el4<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-4/")
el4.seu<-CreateSeuratObject(el4)
el4.seu<-RenameCells(el4.seu, new.names = paste0(unlist(lapply(str_split(colnames(el4.seu), pattern = "-"), "[[", 1)), -4))

el5<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-5/")
el5.seu<-CreateSeuratObject(el5)
el5.seu<-RenameCells(el5.seu, new.names = paste0(unlist(lapply(str_split(colnames(el5.seu), pattern = "-"), "[[", 1)), -5))

el6<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-6/")
el6.seu<-CreateSeuratObject(el6)
el6.seu<-RenameCells(el6.seu, new.names = paste0(unlist(lapply(str_split(colnames(el6.seu), pattern = "-"), "[[", 1)), -6))

el7<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-7/")
el7.seu<-CreateSeuratObject(el7)
el7.seu<-RenameCells(el7.seu, new.names = paste0(unlist(lapply(str_split(colnames(el7.seu), pattern = "-"), "[[", 1)), -7))

el8<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-8/")
el8.seu<-CreateSeuratObject(el8)
el8.seu<-RenameCells(el8.seu, new.names = paste0(unlist(lapply(str_split(colnames(el8.seu), pattern = "-"), "[[", 1)), -8))

el9<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-9/")
el9.seu<-CreateSeuratObject(el9)
el9.seu<-RenameCells(el9.seu, new.names = paste0(unlist(lapply(str_split(colnames(el9.seu), pattern = "-"), "[[", 1)), -9))

el10<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-10/")
el10.seu<-CreateSeuratObject(el10)
el10.seu<-RenameCells(el10.seu, new.names = paste0(unlist(lapply(str_split(colnames(el10.seu), pattern = "-"), "[[", 1)), -10))

el11<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-11/")
el11.seu<-CreateSeuratObject(el11)
el11.seu<-RenameCells(el11.seu, new.names = paste0(unlist(lapply(str_split(colnames(el11.seu), pattern = "-"), "[[", 1)), -11))

el12<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-12/")
el12.seu<-CreateSeuratObject(el12)
el12.seu<-RenameCells(el12.seu, new.names = paste0(unlist(lapply(str_split(colnames(el12.seu), pattern = "-"), "[[", 1)), -12))

el13<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-13/")
el13.seu<-CreateSeuratObject(el13)
el13.seu<-RenameCells(el13.seu, new.names = paste0(unlist(lapply(str_split(colnames(el13.seu), pattern = "-"), "[[", 1)), -13))

el14<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-14/")
el14.seu<-CreateSeuratObject(el14)
el14.seu<-RenameCells(el14.seu, new.names = paste0(unlist(lapply(str_split(colnames(el14.seu), pattern = "-"), "[[", 1)), -14))

el15<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-15/")
el15.seu<-CreateSeuratObject(el15)
el15.seu<-RenameCells(el15.seu, new.names = paste0(unlist(lapply(str_split(colnames(el15.seu), pattern = "-"), "[[", 1)), -15))

el16<-Read10X("/projects/b1042/LinLab/realign-mismatch/output/SoupX/EL-16/")
el16.seu<-CreateSeuratObject(el16)
el16.seu<-RenameCells(el16.seu, new.names = paste0(unlist(lapply(str_split(colnames(el16.seu), pattern = "-"), "[[", 1)), -16))

el.m<-merge(el1.seu, y = c(el2.seu, el3.seu, el4.seu, el5.seu, el6.seu, el7.seu, el8.seu, el9.seu, el10.seu, el11.seu, el12.seu, el13.seu, el14.seu, el15.seu, el16.seu))

#ncol(el.m)
#colnames(el.m)
#colnames(aggr.seu)

#el.m.sorted<-el.m[,colnames(aggr.seu)]
#all(colnames(el.m.sorted) %in% colnames(aggr.seu))
#all(colnames(el.m.sorted) %in% colnames(aggr.seu))

#GetAssayData(aggr.seu, slot = "counts", assay = "RNA")[,"AAACAGCCAAACTGCC-13"]

el.m.sorted<-GetAssayData(el.m, slot = "counts", assay = "RNA")[,colnames(aggr.seu)]
all(colnames(el.m.sorted)==colnames(aggr.seu))

all(rownames(aggr.seu@assays$RNA@counts)==rownames(el.m.sorted))
all(colnames(el.m.sorted)==colnames(aggr.seu))

aggr.seu.replaced<-aggr.seu
aggr.seu.replaced@assays$RNA@counts<-el.m.sorted
aggr.seu.replaced@assays$RNA@data<-el.m.sorted

saveRDS(aggr.seu.replaced, "/projects/b1042/LinLab/realign-mismatch/output/aggr/aggregated-seurat-SoupX.rds")

rm(aggr.seu)
gc()

### Run DoubletFinder
# DoubletFinder.R

### Filter low quality cells using same DoubletFinder QC and remove doublets
aggr.seu.replaced<-readRDS("/projects/b1042/LinLab/realign-mismatch/output/aggr/aggregated-seurat-SoupX.rds")

### No ribo filter
### Remove high nCount_RNA cells (>20000)
aggr.filt2<-subset(aggr.seu.replaced, 
                   subset=percent.mt<1 & 
                     #percent.ribo<5 & 
                     nCount_RNA>350 & 
                     nCount_RNA<20000 &
                     nFeature_RNA > 200 &
                     #nFeature_RNA < 4000 & 
                     nCount_ATAC > 1000 &
                     nCount_ATAC < 100000 &
                     nucleosome_signal < 2 &
                     TSS.enrichment > 2)

### Remove doublets
### Run DoubletFinder.R
doub.ann <- read.table("/projects/b1042/LinLab/realign-mismatch/output/doublets-5percent-relaxed-noRiboFilter-20000nCountRNA.txt", sep="\t", header=T, stringsAsFactors = F, quote="")

aggr.filt2$doublet.annotation<-doub.ann$doublet.ann[match(colnames(aggr.filt2), doub.ann$sample)]

aggr.filt2<-subset(aggr.filt2, subset=doublet.annotation=="Singlet")

saveRDS(aggr.filt2, "/projects/b1042/LinLab/realign-mismatch/output/aggr/aggregated-seurat-SoupX-20000nCountRNA-dedoub.rds")
