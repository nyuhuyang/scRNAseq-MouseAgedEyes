########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gplots)
library(MAST)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file = "data/mm10_young_aged_eyes_2_20190712.Rda"))
Idents(object) <- "Doublets"
object %<>% subset(idents = "Singlet")

object$conditions %<>% factor(levels = c("Young","Aged"))

# Isolate EC
Idents(object) <- "cell.type"
EC <- subset(object, idents = "Endothelial cells")
#========================================
# recluster, Performing CCA integration
#========================================
#======1.4 =========================
set.seed(100)
npcs =75
Idents(EC) = "conditions"
EC_list <- lapply(c("Young","Aged"),function(x) subset(EC,idents=x))
anchors <- FindIntegrationAnchors(EC_list, dims = 1:npcs)
EC <- IntegrateData(anchors, dims = 1:npcs)
remove(anchors,EC_list);GC()
DefaultAssay(EC) <- "integrated"
EC %<>% ScaleData(verbose = FALSE)
EC %<>% RunPCA(npcs = npcs, features = VariableFeatures(EC),verbose = FALSE)
EC %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
EC %<>% FindClusters(reduction = "pca",resolution = 0.3,
                         dims.use = 1:npcs,print.output = FALSE)

meta.data = EC@meta.data[,(ncol(EC@meta.data)-1):ncol(EC@meta.data)]
EC <- subset(object, idents = "Endothelial cells")
EC@meta.data = cbind(EC@meta.data,meta.data[match(rownames(EC@meta.data),rownames(meta.data)),])
EC@meta.data = cbind(EC@meta.data,EC@reductions$umap@cell.embeddings)
EC <- EC[,EC$UMAP_1 <10]
TSNEPlot.1(EC, group.by="integrated_snn_res.0.3",pt.size = 1,label = F,
           label.size = 4, repel = T,title = "All EC cluster in tSNE plot",
           do.print = T)

UMAPPlot.1(EC, group.by="integrated_snn_res.0.3",pt.size = 1,label = F,
           label.size = 4, repel = T,title = "All EC cluster in UMAP plot",
           do.print = T)

#========================================
# Find markers in sub-clusters
#========================================
Idents(EC) = "integrated_snn_res.0.3"
gde.markers <- FindAllMarkers.UMI(EC,logfc.threshold = 0.1,only.pos = FALSE, 
                                       min.pct = 0.1,return.thresh = 0.05)

write.csv(gde.markers,paste0(path,"EC_subcluster_markers.csv"))
(mito.genes <- grep(pattern = "^mt-", x = gde.markers$gene))
if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
GC()
#DoHeatmap.1======
Top_n = 10
top <-  gde.markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)
EC %<>% ScaleData(features= unique(top$gene))

DoHeatmap.1(EC, add.genes = top$gene, Top_n = Top_n, do.print=T, angle = 0,
            group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            label=T, cex.row=7, legend.size = NULL,width=10, height=7,unique.name = T,
            title = paste("Top",Top_n,"DE genes in Endothelial cells sub clusters"))
#========================================
# Find markers between young and age in each cluster
#========================================
EC$cluster_conditions <- paste(EC$conditions,EC$integrated_snn_res.0.3,sep = "_")
Idents(EC) <-"cluster_conditions"

(cell_type <- sort((unique(EC$integrated_snn_res.0.3))))
ident.1 = paste0("Young_",1:5)
ident.2 = paste0("Aged_",1:5)

gde.markers <- FindPairMarkers(EC, ident.1 = c(ident.1,ident.2),
                               ident.2 = c(ident.2,ident.1), only.pos = T,
                               logfc.threshold = 0.05,min.cells.group =3,
                               return.thresh = 0.05,save.files = FALSE)
write.csv(gde.markers,paste0(path,"EC_age_markers.csv"))
(mito.genes <- grep(pattern = "^mt-", x = gde.markers$gene))
if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
GC()
#VlnPlot======
(gde.markers[order(gde.markers$p_val_adj),"gene"][1:10])
EC@meta.data$conditions %<>% factor(levels = c("Young","Aged"))

plots <- VlnPlot(EC, features = c("Gpx3","Ihh"), split.by = "conditions", group.by = "integrated_snn_res.0.3", 
                 pt.size = 0, combine = FALSE)
jpeg(paste0(path,"VlnPlot_EC.jpeg"), units="in", width=10, height=7,res=600)
CombinePlots(plots = plots, ncol = 1)
dev.off()
# Volcano plot=========
(clusters <- unique(gde.markers$cluster1.vs.cluster2))
for(i in 1:(length(clusters)/2)){
        df_young <- gde.markers[gde.markers$cluster1.vs.cluster2 %in% clusters[i],]
        df_aged <- gde.markers[gde.markers$cluster1.vs.cluster2 %in% clusters[i+length(clusters)/2],]

        df_young$log10_p_val_adj = -log10(df_young$p_val_adj)
        df_aged$log10_p_val_adj = -log10(df_aged$p_val_adj)
        df_young$avg_logFC = -df_young$avg_logFC
        df = rbind.data.frame(df_young, df_aged)
        
        df$log10_p_val_adj[df$log10_p_val_adj == "Inf"] = 400
        
        left = rownames(df_young)[df_young$log10_p_val_adj >= head(sort(df_young$log10_p_val_adj,decreasing = T),15) %>% tail(1)]
        right = rownames(df_aged)[df_aged$log10_p_val_adj >= head(sort(df_aged$log10_p_val_adj,decreasing = T),15) %>% tail(1)]
        g <- ggplot(df,aes(avg_logFC,log10_p_val_adj)) + 
                geom_point() + 
                ggtitle(paste("Young (left) and Aged (right) mouse eyes EC's cluster",i)) + 
                ylab("-log10(p_value_adj)")+
                theme_minimal()+
                theme(plot.title = element_text(size=20,hjust = 0.5))+
                ggrepel::geom_text_repel(aes(label = gene), 
                                         data=df[c(left[1:15],right[1:15]),]) +
                geom_point(color = ifelse((df$avg_logFC > 0.05  & df$p_val_adj < 0.05) , "red",
                                          ifelse((df$avg_logFC < -0.05 & df$p_val_adj < 0.05), "blue","gray")))
        jpeg(paste0(path,"Volcano_plot_EC_",clusters[i],".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
}

