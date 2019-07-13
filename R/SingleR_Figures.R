library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
source("R/utils/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/mm10_young_aged_eyes_2_20190712.Rda"))
(load(file = "output/singler_T_mm10_young_aged_eyes_2_20190712.Rda"))
# if singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(object@assays$RNA@data)
if(length(singler$singler[[1]]$SingleR.single$labels) < ncol(object@assays$RNA@data)){
        all.cell = colnames(object);length(all.cell)
        know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
        object = object[,know.cell]
}


table(rownames(singler$singler[[1]]$SingleR.single$labels) %in% colnames(object))
singler$meta.data$orig.ident = object@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = object@reductions$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = Idents(object) # the Seurat clusters (if 'clusters' not provided)
save(singler,file="output/singler_T_mm10_young_aged_eyes_2_20190712.Rda")

##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                       "singler1main" = singler$singler[[1]]$SingleR.single.main$labels,
                       "orig.ident" = object@meta.data$orig.ident,
                       "res.0.6" = object@meta.data$integrated_snn_res.0.6,
                       row.names = colnames(object))

table(rownames(singlerDF) %in% colnames(object))
head(singlerDF)
apply(singlerDF,2,function(x) length(unique(x)))
##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,normalize = F))
dev.off()
jpeg(paste0(path,"DrawHeatmap_sub1_N.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single,top.n = 50,normalize = T))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singlerDF$singler1sub, singlerDF$orig.ident)) %>%
        kable_styling()
singlerDF$orig.ident %>% table() %>% kable() %>% kable_styling()
singlerDF$cell.type %>% table() %>% kable() %>% kable_styling()

singlerDF$singler1sub = gsub("B cells|T cells|NK cells",
                             "Lymphoid cells",singlerDF$singler1sub)
singlerDF$singler1sub = gsub("Dendritic cells|Granulocytes|Macrophages|Macrophages activated|Monocytes",
                             "Myeloid cells",singlerDF$singler1sub)
singlerDF$singler1sub = gsub("Fibroblasts|Fibroblasts activated|Fibroblasts senescent",
                             "Fibroblasts",singlerDF$singler1sub)
singlerDF$cell.type <- plyr::mapvalues(singlerDF$res.0.6,
                                         from = 0:18,
                                           to = c("Stromal cells",
                                                  "Endothelial cells",
                                                  "Retinal pigment epithelium",
                                                  "Retinal pigment epithelium",
                                                  "Pericytes",
                                                  "Pericytes",
                                                  "Endothelial cells",
                                                  "Endothelial cells",
                                                  "Stromal cells",
                                                  "Myeloid cells",
                                                  "Retinal pigment epithelium",
                                                  "Endothelial cells",
                                                  "Myelinating Schwann cells",
                                                  "Myeloid cells",
                                                  "Lymphoid cells",
                                                  "Stromal cells",
                                                  "Pericytes",
                                                  "Myelinating Schwann cells",
                                                  "Melanocytes"))
##############################
# process color scheme
##############################

singler_colors <- readxl::read_excel("doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors1[duplicated(singler_colors1)]
length(singler_colors1)
apply(singlerDF[,c("singler1sub","cell.type")],2,function(x) length(unique(x)))
singlerDF[,c("cell.type")] %>% table() %>% kable() %>% kable_styling()
object <- AddMetaData(object = object,metadata = singlerDF)
object <- AddMetaColor(object = object, label= "cell.type", colors = singler_colors1)
Idents(object) <- "cell.type"
object %<>% sortIdent()
# Idents(object) <- "Doublets"
#object %<>% subset(idents = "Singlet")
TSNEPlot.1(object, cols = ExtractMetaColor(object),label = F,pt.size = 1,no.legend = T,
         label.size = 5, repel = T,do.print = T,title = "All cell types in tSNE plot")

UMAPPlot.1(object, group.by="cell.type",pt.size = 1,label = F,
         cols = ExtractMetaColor(object),do.print = T,no.legend = T,
         label.size = 5, repel = T,title = "All cell types in UMAP plot")

save(object,file="data/mm10_young_aged_eyes_2_20190712.Rda")

