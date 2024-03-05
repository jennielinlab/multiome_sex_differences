### Call MACS2

library(tidyverse)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(GenomicFeatures)
library(EnsDb.Mmusculus.v79)

#annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
#seqlevelsStyle(annotation) <- "UCSC"
annotation<-readRDS("/projects/b1042/LinLab/realign-mismatch/data/EnsDb.Mmusculus.v79-annotation.rds")

# <- readRDS("/projects/b1042/LinLab/realign-mismatch/output/aggr/aggr-filtered-and-integrated-harmony.rds")
aggr.harm<-readRDS("/projects/b1042/LinLab/realign-mismatch/output/aggr/aggregated-seurat-SoupX-20000nCountRNA-dedoub-harmony.rds")

fragpath<-"/projects/b1042/LinLab/realign-mismatch/output/aggr/NTS_combined/outs/atac_fragments.tsv.gz"

# call peaks using MACS2
Idents(aggr.harm)<-"cellType"
DefaultAssay(aggr.harm)<-"ATAC"

peaks <- CallPeaks(aggr.harm, macs2.path = "/software/MACS2/bin/macs2", group.by = "cellType")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(aggr.harm),
  features = peaks,
  cells = colnames(aggr.harm)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
aggr.harm[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation,
  sep = c(":", "-"),
  genome = 'mm10'
)

DefaultAssay(aggr.harm)<-"peaks"

### Remove ATAC assay
aggr.harm[["ATAC"]]<-NULL

saveRDS(aggr.harm, "/projects/b1042/LinLab/realign-mismatch/output/aggr/aggregated-seurat-SoupX-20000nCountRNA-dedoub-harmony-MACS-noATAC-celltype.rds")
