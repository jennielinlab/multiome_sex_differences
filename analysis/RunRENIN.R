### Run RENIN
### Matt Wright, adapted from https://nledru.github.io/RENIN/gettingstarted.html
### last updated: October 4, 2023

library(tidyverse)
library(Seurat)
library(Signac)
library(RENIN)
library(VISION)
library(chromVARmotifs)
library(BSgenome.Mmusculus.UCSC.mm10)
library(harmony)

### Define cell types of interest
coi<-c("PTS1", "PTS2", "PTS3", "Endothelial", "Fibroblast", "LOH", "Podocytes")

#aggr.filt<-readRDS("/media/data/projects/Matt/NTS-sex-differences/output/aggregated-20000nCountRNA-noMT-genes-noEndoPTcluster.rds")

### Add new mouse motifs
#data("mouse_pwms_v2") # filtered human motifs from cisBP
#motifs <- mouse_pwms_v2
#aggr.filt <- AddMotifs(aggr.filt, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = motifs, assay = "peaks")
#saveRDS(aggr.filt, "/media/data/projects/Matt/NTS-sex-differences/output/aggregated-20000nCountRNA-noMT-genes-noEndoPTcluster-with-RENIN-motifs.rds")

### Read in Seurat
aggr.filt<-readRDS("/media/data/projects/Matt/NTS-sex-differences/output/aggregated-20000nCountRNA-noMT-genes-noEndoPTcluster-with-RENIN-motifs.rds")

DefaultAssay(aggr.filt)<-"RNA"
aggr.filt<-NormalizeData(aggr.filt)

allgenes<-rownames(aggr.filt)
MTgenes<-allgenes[grepl(allgenes, pattern="^mt-")]
noMTgenes<-allgenes[!grepl(allgenes, pattern="^mt-")]

DefaultAssay(aggr.filt)<-"peaks"

Idents(aggr.filt)<-"cellType2"

#aggr.filt<-subset(aggr.filt, features=noMTgenes)

### 
mats <- prepare_pseudocell_matrix(aggr.filt, assay = c("peaks", "RNA"), cells_per_partition = 100, reduction1 = "harmony.rna", reduction2 = "harmony.peaks")
expr_mat <- mats[["RNA"]]
peak_mat <- mats[["peaks"]]

### Remove MT genes from expression matrix
#expr_mat<-expr_mat[,-which(colnames(expr_mat) %in% MTgenes)]


comparisons<-c("MvF", "YvA", "CvN")
#comparisons<-c("YvA", "CvN")
for(comp in comparisons){
  for(i in coi){
    celltype<-i

    if(comp=="MvF"){
      c1<-c("Young", "Aged")
      c2<-c("Control", "NTS")
    } 
    if(comp=="CvN") {
      c1<-c("M", "F")
      c2<-c("Young", "Aged")
    } 
    if(comp=="YvA"){
      c1<-c("M", "F")
      c2<-c("Control", "NTS")
    }
    
    for(j in c1){
      for(k in c2){
        
        cat(paste0("Running: ", comp, ", ", j, ", ", k, ", in ", celltype, "\n"))
        
        if(celltype=="Endothelial" & comp=="CvN" & j=="M" & k=="Aged"){
          # thank u
          next
        }
          
        if(comp=="MvF"){
          id1<-paste0("F.", j, ".", k)
          id2<-paste0("M.", j, ".", k)
          up="Female"
          down="Male"
        }
        if(comp=="CvN"){
          id1<-paste0(j, ".", k, ".Control")
          id2<-paste0(j, ".", k, ".NTS")
          up="Control"
          down="NTS"
        }
        if(comp=="YvA"){
          id1<-paste0(j, ".Young.", k)
          id2<-paste0(j, ".Aged.", k)
          up="Young"
          down="Aged"
        }
    
        ### MvF, Young, Control
        ct.sub<-subset(aggr.filt, subset=cellType %in% celltype )
        ct.sub<-subset(ct.sub, subset=sex.age.treat %in% c(id1, id2))
      
        de.genes<-FindMarkers(aggr.filt, ident.1 = id1, ident.2 = id2, assay = "RNA", subset.ident = celltype, group.by = "sex.age.treat", test.use = "LR", features = noMTgenes)
        #de.genes.sig<-de.genes[which(de.genes$p_val_adj<0.05),]
        gene_list<-rownames(de.genes)
        
        ### Identify putative cis-regulatory elements for each of the input genes (DE genes)
        peak_results <- run_peak_aen(aggr.filt, expr_mat, peak_mat, gene_list, lambda2 = 0.5, max_distance = 5e+05, num_bootstraps = 100)
        aen_lists <- make_aen_lists(peak_results)
        
        
        up_genes <- rownames(de.genes)[which(de.genes$avg_log2FC > 0)]
        down_genes <- rownames(de.genes)[which(de.genes$avg_log2FC < 0)]
        all_cres <- unique(unlist(aen_lists))
        
        
        cre_scores <- lapply(peak_results, function(x) x[[4]][union(1, which(x[[4]][, "coef_if_kept"] != 0)), "coef_if_kept"] * ifelse(x[[1]] %in% down_genes, -1, 1))
        cre_scores <- cre_scores[which(unlist(lapply(cre_scores, length)) > 1)]
        
        cre_scores.out<-list()
        
        for(i in 1:length(cre_scores)){
          genename=names(cre_scores)[i]
          
          tmpcoefs<-as.data.frame(cre_scores[[i]])
          colnames(tmpcoefs)<-"CRE.Score"
          
          tmpcoefs %>% rownames_to_column("Peak") -> tmpcoefs
          tmpcoefs$DEgene<-genename
          
          tmpcoefs<-tmpcoefs[which(tmpcoefs$CRE.Score!=0),]
          tmpcoefs<-tmpcoefs[which(tmpcoefs$Peak!="(Intercept)"),]
          
          cre_scores.out[[genename]]<-tmpcoefs
        }
        
        cre_scores.out.df<-bind_rows(cre_scores.out)
        
        #cre_scores.out.df$condition<-ifelse(cre_scores.out.df$CRE.Score>0, up, down)
        #cre_scores.out.df$condition<-ifelse(cre_scores.out.df$CRE.Score==0, "NS", cre_scores.out.df$condition)
        
        outname<-paste0("/media/data/projects/Matt/NTS-sex-differences/output/RENIN/noMT/RNA/allGenes/", comp, "/", j, "/", k, "/Peak-gene-results/", celltype, ".txt")
        write.table(cre_scores.out.df, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
        
        cre_total_scores <- bind_rows(cre_scores)
        cre_total_scores[is.na(cre_total_scores)] <- 0
        cre_total_scores <- cre_total_scores[,-1]
        cre_total_scores_sums <- colSums(cre_total_scores)
        cre_num_genes <- apply(cre_total_scores, 2, function(x) length(which(x != 0))) # how many genes each CRE targets
        up_cres <- names(cre_total_scores_sums)[which(cre_total_scores_sums > 0)]
        down_cres <- names(cre_total_scores_sums)[which(cre_total_scores_sums < 0)]
        
        #pdf("../plots/DE/20000nCountRNA/RENIN/part1-M-TFs-test-v2.pdf", width = 8, height = 6)
        #plot_motif_enrichment(aggr.filt, m_cres, num_top_label = 20)
        #dev.off()
        
        m.cre.motifs<-FindMotifs(object = aggr.filt, features = up_cres)
        outname<-paste0("/media/data/projects/Matt/NTS-sex-differences/output/RENIN/noMT/RNA/allGenes/", comp, "/", j, "/", k, "/Motifs/", celltype, "-", up, "-Motifs.txt")
        write.table(m.cre.motifs, file = outname, sep = "\t", col.names = T, row.names = F, quote=F)
        
        
        f.cre.motifs<-FindMotifs(object = aggr.filt, features = down_cres)
        outname<-paste0("/media/data/projects/Matt/NTS-sex-differences/output/RENIN/noMT/RNA/allGenes/", comp, "/", j, "/", k, "/Motifs/", celltype, "-", down, "-Motifs.txt")
        write.table(f.cre.motifs, file = outname, sep = "\t", col.names = T, row.names = F, quote=F)
      
        ### PART 2
        expr_mat.ct <- prepare_pseudocell_matrix(ct.sub, assay = "RNA", cells_per_partition = 10, reduction1="harmony.rna", reduction2="harmony.peaks")
        ### Remove MT genes from expression matrix
        expr_mat.ct<-expr_mat.ct[,-which(colnames(expr_mat.ct) %in% MTgenes)]
        
        
        tf_results <- run_tf_aen(ct.sub, expr_mat.ct, peak_results, gene_list, lambda2 = 0.5)
        
        regulator_tf_names <- unlist(aggr.filt@assays$peaks@motifs@motif.names)
        regulator_tf_names <- regulator_tf_names[which(regulator_tf_names %in% rownames(GetAssayData(aggr.filt, assay = "RNA")))]
        
        tf.results.out<-list()
        
        for(i in 1:length(tf_results)){
          genename=names(tf_results)[i]
          
          tmpcoefs<-as.data.frame(tf_results[[i]]$coefs)
          
          tmpcoefs %>% rownames_to_column("TF") -> tmpcoefs
          tmpcoefs$DEgene<-genename
          
          tmpcoefs<-tmpcoefs[which(tmpcoefs$coef_if_kept!=0),]
          
          tf.results.out[[genename]]<-tmpcoefs
        }
        
        tf.results.out.df<-bind_rows(tf.results.out)
        
        tf.results.out.df<-tf.results.out.df[,c(5,1,2,4)]
        colnames(tf.results.out.df)<-c("DEgene", "RegulatoryTF", "RegulatoryWeight", "SE")
        
        tf.results.out.df<-tf.results.out.df[which(tf.results.out.df$RegulatoryTF %in% regulator_tf_names),]
        #tf.results.out.df$condition<-ifelse(tf.results.out.df$coef_if_kept>0, up, down)
        #tf.results.out.df$condition<-ifelse(tf.results.out.df$coef_if_kept==0, "NS", tf.results.out.df$condition)
        
        outname<-paste0("/media/data/projects/Matt/NTS-sex-differences/output/RENIN/noMT/RNA/allGenes/", comp, "/", j, "/", k, "/TF-gene-results/", celltype, ".txt")
        write.table(tf.results.out.df, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
        
        
        ### Rank transcription factor activities across all genes. Weight downregulated genes negatively
        tf_df <- rank_tfs(tf_results, ct.sub, down_genes, regulator_tf_names, num_cores = 1)
        tf_df$condition<-ifelse(tf_df$Score>0, up, down)
        tf_df$condition<-ifelse(tf_df$Score==0, "NS", tf_df$condition)
        
        outname<-paste0("/media/data/projects/Matt/NTS-sex-differences/output/RENIN/noMT/RNA/allGenes/", comp, "/", j, "/", k, "/", celltype, ".txt")
        write.table(tf_df, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
            
        #tf_df
        #head(tf_df)
        
        #pdf("../plots/DE/20000nCountRNA/RENIN/part2-test-v2.pdf", width = 8, height = 6)
        #plot_tf_rankings(tf_df)
        #+
        #  scale_fill_discrete(labels=c("F", "M"))+
        #  labs(fill="Sex")
        #dev.off()
        
        centrality_rankings <- rank_tfs_by_centrality(tf_results, ct.sub)
        
        centrality_rankings.out<-as.data.frame(centrality_rankings$Betweenness)
        colnames(centrality_rankings.out)<-"Betweenness"
        centrality_rankings.out %>%
           rownames_to_column("TF") -> centrality_rankings.out
        outname<-paste0("/media/data/projects/Matt/NTS-sex-differences/output/RENIN/noMT/RNA/allGenes/", comp, "/", j, "/", k, "/", celltype, "-Betweenness.txt")
        write.table(centrality_rankings.out, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
        
        centrality_rankings.out<-as.data.frame(centrality_rankings$PageRank)
        colnames(centrality_rankings.out)<-"PageRank"
        centrality_rankings.out %>%
          rownames_to_column("TF") -> centrality_rankings.out
        outname<-paste0("/media/data/projects/Matt/NTS-sex-differences/output/RENIN/noMT/RNA/allGenes/", comp, "/", j, "/", k, "/", celltype, "-PageRank.txt")
        write.table(centrality_rankings.out, file = outname, sep="\t", col.names = T, row.names = F, quote=F)
        
        #head(centrality_rankings[[1]], n = 10)
        #head(centrality_rankings[[2]], n = 10)
        
        rm(ct.sub)
        gc()
    
      }
    }
  }
}
