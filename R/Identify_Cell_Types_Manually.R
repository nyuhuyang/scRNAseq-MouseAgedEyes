library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
source("R/utils/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#====== 2.1 pathway analysis ==========================================
(load(file="data/mm10_young_aged_eyes_2_20190712.Rda"))
DefaultAssay(object) <- "RNA"
df_markers <- readxl::read_excel("doc/Celltype Markers.xlsx",sheet = "Sheet1")

#df_markers <- readxl::read_excel("../seurat_resources/bio-rad-markers.xlsx",sheet = "Human.sub")
colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
markers = df_markers[,-grep("Alias",colnames(df_markers))]
marker.list <- df2list(markers)

marker.list %<>% lapply(function(x) x[1:16]) %>% 
     lapply(function(x) FilterGenes(object,x)) %>% 
     lapply(function(x) x[!is.na(x)]) %>% 
    lapply(function(x) x[1:min(length(x),9)])
marker.list %>% list2df %>% t %>% kable() %>% kable_styling()

dir.create(paste0(path,"markers"), recursive = T)
Idents(object) = "integrated_snn_res.0.6"
for(i in 1:length(marker.list)){
    if(length(marker.list[[i]]) == 0) next
    p <- lapply(marker.list[[i]], function(marker) {
        FeaturePlot(object = object, feature = marker,pt.size = 0.5,
                    reduction="tsne", label = T)+
            NoLegend()+
            ggtitle(paste0(marker,Alias(df = df_markers,gene = marker)))+
            theme(plot.title = element_text(hjust = 0.5,size = 15,face = "plain"))
    })
    jpeg(paste0(path,"markers/",names(marker.list)[i],".jpeg"),units="in", width=10, height=7,res=600)
    print(do.call(cowplot::plot_grid, p)+ ggtitle(paste(names(marker.list)[i],"markers"))+
              theme(plot.title = element_text(hjust = 0.5,size = 20)))
    dev.off()
    print(paste0(i,":",length(marker.list)))
}

object_data <- object@assays$RNA@data
save(object_data, file = "data/object_mm10_young_aged_eyes_2_20190712.Rda")

#============================
# Dot plot
#============================
markers.to.plot <- c("Cdh5", "Pecam1", "Dcn", "Pdgfra", "Acta2", "Myh11", 
                     "Pdgfrb","Rgs5","Rpe65","Rlbp1","Pmel","Mlana","Mbp","Mpz",
                     "Ptprc", "Laptm5", "Cd14", "Ms4a7","Nkg7","Cd3g","Cd19")
markers.to.plot <- FilterGenes(object,markers.to.plot)
Idents(object) = "res.0.6"
Idents(object) %<>% factor(levels = c(9,13,14,12,17,18,2,3,10,4,5,16,0,8,15,1,6,7,11))
g <- DotPlot(object, features = rev(markers.to.plot), 
             cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
jpeg(paste0(path,"Dot_plot_cluster.jpeg"), units="in", width=10, height=6,res=600)
print(g)
dev.off()


Idents(object) = "cell.type"
Idents(object) %<>% factor(levels = c("Lymphoid cells","Myeloid cells",
                                      "Myelinating Schwann cells","Melanocytes",
                                      "Retinal pigment epithelium","Pericytes","Stromal cells",
                                      "Endothelial cells"))
g <- DotPlot(object, features = rev(markers.to.plot), 
             split.by = "conditions",
             cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
jpeg(paste0(path,"Dot_plot.jpeg"), units="in", width=10, height=7,res=600)
print(g)
dev.off()
