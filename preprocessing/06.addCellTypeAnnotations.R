### Add celltype annotations
aggr.filt<-readRDS("../output/aggr/aggregated-seurat-SoupX-20000nCountRNA-dedoub-harmony.rds")

aggr.filt<-FindSubCluster(aggr.filt, cluster = 19, graph.name = "wsnn", subcluster.name = "Clust19")
aggr.filt<-FindSubCluster(aggr.filt, cluster = 29, graph.name = "wsnn", subcluster.name = "Clust29", resolution = 0.8)
aggr.filt<-FindSubCluster(aggr.filt, cluster = 18, graph.name = "wsnn", subcluster.name = "Clust18", resolution = 0.8)

ct.tbl<-as.data.frame(table(aggr.filt$predictedCellType, aggr.filt$Clust19))
ct.tbl %>%
  pivot_wider(names_from = "Var1", values_from = "Freq") -> ct.tbl
ct.tbl$predictedCT<-colnames(ct.tbl[,2:ncol(ct.tbl)])[apply(ct.tbl[,2:ncol(ct.tbl)], MARGIN = 1, FUN = which.max)]


VlnPlot(aggr.filt, features = "Top2a", group.by = "Clust29")


aggr.filt$cellType<-aggr.filt$predictedCellType


#aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% 19, "Macrophages", aggr.filt$cellType)
#aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% 29, "cyclingPT", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$Clust29 %in% c("29_0", "29_2"), "cyclingPT", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$Clust29 %in% c("29_1", "29_3"), "PTS3", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% 15, "ICB", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% 23, "ICA", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% 13, "Podocytes", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(9, 14, 20), "Endothelial", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(26), "PEC", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(10), "Fibroblast", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(31), "Pericytes", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(8), "InjuredPT", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(28), "MaculaDensa", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(4, 5, 18), "LOH", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(22), "InjuredLOH", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(24), "ThinLimb", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$Clust19 %in% c("19_0", "19_2"), "Macrophages", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$Clust19 %in% c("19_3"), "PTS1", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$Clust19 %in% c("19_1", "19_4"), "PTS2", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(1, 16, 17), "PTS1", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(0), "PTS2", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(2, 6, 25, 30), "PTS3", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(3, 27), "DCT", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$Clust18 %in% c("18_2"), "DCT", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(7), "DCT.CNT", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(12, 21), "CNT", aggr.filt$cellType)
aggr.filt$cellType<-ifelse(aggr.filt$seurat_clusters %in% c(11), "PC", aggr.filt$cellType)

aggr.filt<-FindSubCluster(aggr.filt, cluster = 18, graph.name = "wsnn", subcluster.name = "Clust18", resolution = 0.2)
DimPlot(aggr.filt, group.by = "Clust18", reduction = "wnn.umap", label = T)+NoLegend()

aggr.filt$cellType<-ifelse(aggr.filt$Clust18 %in% c("18_2"), "DCT", aggr.filt$cellType)


DimPlot(aggr.filt, reduction = "wnn.umap", label = T)+NoLegend()
DimPlot(aggr.filt, group.by = "cellType", reduction = "wnn.umap", label = T)+NoLegend()
DimPlot(aggr.filt, group.by = "Phase", reduction = "wnn.umap")
DimPlot(aggr.filt, group.by = "Clust19", reduction = "wnn.umap", label = T)+NoLegend()
DimPlot(aggr.filt, group.by = "predictedCellType", reduction = "wnn.umap", label = T)+NoLegend()

aggr.filt$cellType2<-factor(aggr.filt$cellType, levels = rev(ct.order.names))



DotPlot(aggr.filt, features = c("Lrp2", "Slc5a12", "Slc13a3", "Slc6a18",  "Umod", "Slc12a3", "Slc8a1", "Egfem1", "Cfh", "Flt1", "Aqp2", "Nphs2", "Ncam1", "Ebf1", "Aqp6", "Slc26a4", "Slc14a2", "Epha7", "Cd74", "Nos1", "Thsd4", "Il34", "Vcam1", "Havcr1", "Lcn2", "Top2a"), group.by = "cellType2", scale.max = 80)+
  scale_color_viridis_c(direction = -1)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        plot.background = element_rect(fill = 'white'))
ggsave("../plots/clustering/20230626-DotPlot-v7.pdf", width = 12, height = 8, units = "in")

table(aggr.filt$cellType)


saveRDS(aggr.filt, "../output/aggr/aggregated-seurat-SoupX-20000nCountRNA-dedoub-harmony.rds")
