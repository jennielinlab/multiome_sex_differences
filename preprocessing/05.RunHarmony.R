library(tidyverse)
library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomicFeatures)
library(DoubletFinder)
library(EnsDb.Mmusculus.v79)
library(harmony)

#library(AnnotationHub)
#ah = AnnotationHub()

### Find annotation that matches the cellranger arc build (homo sapiens v98)
#annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
#seqlevelsStyle(annotation) <- "UCSC"
#saveRDS(annotation, "/projects/b1042/LinLab/realign-mismatch/data/EnsDb.Mmusculus.v79-annotation.rds")

annotation<-readRDS("/projects/b1042/LinLab/realign-mismatch/data/EnsDb.Mmusculus.v79-annotation.rds")

chr.y<-annotation[seqnames(annotation)=="chrY"]

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

### Add metadata and run harmony on dedoub sample

aggr.filt <- readRDS("../output/aggr/aggregated-seurat-SoupX-20000nCountRNA-dedoub.rds")

aggr.filt$sex<-ifelse(aggr.filt$sample %in% c(1, 10, 13, 14, 2, 5, 6, 9), "F", "M")
aggr.filt$age<-ifelse(aggr.filt$sample %in% c(1, 3, 5, 7, 9, 11, 13, 15), "Young", "Aged")
aggr.filt$treatment<-ifelse(aggr.filt$sample %in% c(1:8), "Control", "NTS")
#aggr.filt$sex<-ifelse(aggr.filt$sample %in% c("EL-1", "EL-10", "EL-13", "EL-14", "EL-2", "EL-5", "EL-6", "EL-9"), "F", "M")
#aggr.filt$age<-ifelse(aggr.filt$sample %in% c("EL-1", "EL-3", "EL-5", "EL-7", "EL-9", "EL-11", "EL-13", "EL-15"), "Young", "Aged")
#aggr.filt$treatment<-ifelse(aggr.filt$sample %in% paste0("EL-", 1:8), "Control", "NTS")
aggr.filt$sex.age<-paste0(aggr.filt$sex, ".", aggr.filt$age)
aggr.filt$sex.age.treat<-paste0(aggr.filt$sex, ".", aggr.filt$age, ".", aggr.filt$treatment)

allmd<-read.table("../output/aggr/allbarcode-metrics.txt", sep="\t", header=T, stringsAsFactors=F, quote="")

allmd$barcode2<-paste0(unlist(lapply(str_split(allmd$barcode, pattern = "-"), "[[", 1)), "-",
                       unlist(lapply(str_split(allmd$sample, pattern = "-"), "[[", 2)))

#md.filt<-allmd[which(allmd$barcode2 %in% colnames(aggr.harm)),]
#md.filt<-md.filt[match(colnames(aggr.harm), md.filt$barcode2),]
#md.filt<-md.filt[match(md.filt$barcode, colnames(aggr.harm)),]
#all(md.filt$barcode2==colnames(aggr.harm))
allmd$sample<-NULL

rownames(allmd)<-allmd$barcode2

aggr.filt<-AddMetaData(aggr.filt, metadata = allmd)
aggr.filt$barcode2<-aggr.filt$barcode<-NULL


DefaultAssay(aggr.filt)<-"RNA"
aggr.filt<-NormalizeData(aggr.filt)

aggr.filt<-CellCycleScoring(object = aggr.filt, s.features = m.s.genes, g2m.features = m.g2m.genes)

#VlnPlot(aggr.filt, features = "Xist", group.by = "sex")

#aggr.filt<-SCTransform(aggr.filt, vst.flavor = "v2", verbose = F, vars.to.regress = c("percent.mt", "percent.rps", "percent.rpl", "S.Score", "G2M.Score"))
aggr.filt<-SCTransform(aggr.filt, vst.flavor = "v2", verbose = F, vars.to.regress = c("percent.mt", "percent.ribo"))

# These are now standard steps in the Seurat workflow for visualization and clustering
aggr.filt <- RunPCA(aggr.filt, verbose = FALSE)

set.seed(42)
aggr.filt<-RunHarmony(aggr.filt, group.by.vars="sample", reduction="pca", assay.use="SCT", reduction.save="harmony.rna")
aggr.filt<-RunUMAP(aggr.filt, reduction="harmony.rna", dims = 1:30, reduction.name="rna.umap", reduction.key="rnaUMAP_")

DimPlot(aggr.filt, reduction = "rna.umap", group.by = "sex.age.treat")
DimPlot(aggr.filt, reduction = "rna.umap", group.by = "sample")



DefaultAssay(aggr.filt) <- "ATAC"
aggr.filt <- FindTopFeatures(aggr.filt, min.cutoff = 'q5')
aggr.filt <- RunTFIDF(aggr.filt)
aggr.filt <- RunSVD(aggr.filt)

set.seed(42)
aggr.filt<-RunHarmony(aggr.filt, group.by.vars="sample", reduction="lsi", assay.use="ATAC", reduction.save="harmony.atac", project.dim=F)
aggr.filt<-RunUMAP(aggr.filt, reduction="harmony.atac", dims = 2:30, reduction.name="atac.umap", reduction.key="atacUMAP_")



aggr.filt<-FindMultiModalNeighbors(aggr.filt, reduction.list=list("harmony.rna","harmony.atac"), dims.list=list(1:30,2:50))
aggr.filt<-RunUMAP(aggr.filt, nn.name="weighted.nn", reduction.name="wnn.umap", reduction.key="wnnUMAP_")
aggr.filt<-FindClusters(aggr.filt, resolution = 1, graph.name = "wsnn", )
DefaultAssay(aggr.filt)<-"SCT"

DimPlot(aggr.filt, reduction="rna.umap")
DimPlot(aggr.filt, reduction="atac.umap")
DimPlot(aggr.filt, reduction = "wnn.umap")

DimPlot(aggr.filt, reduction = "wnn.umap", group.by = "sex.age.treat")
DimPlot(aggr.filt, reduction = "wnn.umap", group.by = "sample")
DimPlot(aggr.filt, reduction = "wnn.umap", group.by = "sex")
DimPlot(aggr.filt, reduction = "wnn.umap", group.by = "age")
DimPlot(aggr.filt, reduction = "wnn.umap", group.by = "treatment")


DimPlot(aggr.filt, reduction = "wnn.umap", label = T)+NoLegend()
DimPlot(aggr.filt, reduction = "wnn.umap", group.by = "Phase")

saveRDS(aggr.filt, "../output/aggr/aggregated-seurat-SoupX-20000nCountRNA-dedoub-harmony.rds")
